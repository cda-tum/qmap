//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"

#include "Logic.hpp"
#include "ir/operations/OpType.hpp"
#include "sc/utils.hpp"

#include <cstddef>
#include <plog/Log.h>
#include <plog/Severity.h>
#include <utility>

namespace cs::encoding {

using namespace logicbase;

void SingleGateEncoder::assertConsistency() const {
  PLOG_DEBUG << "Asserting gate consistency";
  LogicVector gateVariables{};
  gateVariables.reserve(N * (1 + SINGLE_QUBIT_GATES.size()));
  for (std::size_t t = 0U; t < T; ++t) {
    for (std::size_t q = 0U; q < N; ++q) {
      vars.collectSingleQubitGateVariables(t, q, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, true, gateVariables);
    }
    IF_PLOG(plog::verbose) {
      PLOG_VERBOSE << "Gate variables at time " << t;
      for (const auto& var : gateVariables) {
        PLOG_VERBOSE << var.getName();
      }
    }
    assertExactlyOne(gateVariables);
    gateVariables.clear();
  }
}

void SingleGateEncoder::assertGateConstraints() {
  PLOG_DEBUG << "Asserting gate constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    PLOG_VERBOSE << "Asserting gate constraints at time " << t;
    assertSingleQubitGateConstraints(t);
    assertTwoQubitGateConstraints(t);
    assertNoGateNoChangeConstraint(t);
  }
}

void SingleGateEncoder::assertNoGateNoChangeConstraint(const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    const auto noChange = createNoChangeOnQubit(pos, q);
    const auto noGate = createNoGateOnQubit(pos, q);
    lb->assertFormula(LogicTerm::implies(noGate, noChange));
  }
}

void SingleGateEncoder::assertSingleQubitGateConstraints(std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    PLOG_DEBUG << "Asserting gates on " << q;
    assertZConstraints(pos, q);
    assertXConstraints(pos, q);
    assertRConstraints(pos, q);
  }
}

void SingleGateEncoder::assertTwoQubitGateConstraints(const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      if (couplingMap.find(Edge{ctrl, trgt}) == couplingMap.end()) {
        PLOG_DEBUG << "Asserting no CNOT on " << ctrl << " and " << trgt;
        lb->assertFormula(!twoQubitGates[ctrl][trgt]);
        continue;
      }
      const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);

      PLOG_DEBUG << "Asserting CNOT on " << ctrl << " and " << trgt;

      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
    }
  }
}

LogicTerm SingleGateEncoder::createTwoQubitGateConstraint(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  auto changes = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos + 1, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos + 1, ctrl, trgt);

  changes = changes && (tvars->x[pos + 2][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 2][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 2][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 2][trgt] == zTrgt);
  changes = changes &&
            (tvars->r[pos + 2] ==
             (tvars->r[pos + 1] ^ tvars->twoQubitRChange(pos + 1, ctrl, trgt)));

  return changes;
}

LogicTerm SingleGateEncoder::createNoChangeOnQubit(const std::size_t pos,
                                                   const std::size_t q) {
  auto noChange = LogicTerm(true);
  noChange = noChange && (tvars->x[pos + 2][q] == tvars->x[pos + 1][q]);
  noChange = noChange && (tvars->z[pos + 2][q] == tvars->z[pos + 1][q]);
  return noChange;
}

LogicTerm SingleGateEncoder::createNoGateOnQubit(const std::size_t pos,
                                                 const std::size_t q) {
  const auto& singleQubitGates = vars.gS[pos];
  auto noGate = LogicTerm(true);
  for (const auto& gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    noGate = noGate && !singleQubitGates[gateToIndex(gate)][q];
  }
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t i = 0; i < N; ++i) {
    if (i == q) {
      continue;
    }
    noGate = noGate && !twoQubitGates[i][q];
    noGate = noGate && !twoQubitGates[q][i];
  }

  return noGate;
}

void SingleGateEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1];

  // collect variables of single-qubit gates that could be applied to `qubit`
  auto singleQubitGate = LogicTerm(false);
  for (const auto& gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    // any single-qubit gate on qubit q
    singleQubitGate = singleQubitGate || gSNow[gateToIndex(gate)][qubit];
  }

  // collect gate variables of the next timestep that should not be applied.
  auto disallowed = LogicTerm(true);

  // no single-qubit gate on a lower qubit
  for (std::size_t lower = 0U; lower < qubit; ++lower) {
    for (const auto& gate : SINGLE_QUBIT_GATES) {
      if (gate == qc::OpType::None) {
        continue;
      }
      disallowed = disallowed && !gSNext[gateToIndex(gate)][lower];
    }
  }

  lb->assertFormula(LogicTerm::implies(singleQubitGate, disallowed));

  // once no gate is applied, no other gate can be applied, i.e., in the next
  // timestep any of the `None` gate variables must be selected.
  const auto noneIndex = gateToIndex(qc::OpType::None);
  auto noGate = LogicTerm(false);
  for (std::size_t q = 0U; q < N; ++q) {
    noGate = noGate || gSNext[noneIndex][q];
  }
  lb->assertFormula(LogicTerm::implies(gSNow[noneIndex][qubit], noGate));
}

void SingleGateEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gCNext = vars.gC[pos + 1];
  const auto& gSNext = vars.gS[pos + 1];
  for (const auto& [control, target] :
       {std::pair{ctrl, trgt}, std::pair{trgt, ctrl}}) {
    const auto& current = vars.gC[pos][control][target];

    // two identical CNOTs may not be applied in a row because they would
    // cancel.
    auto disallowed = !gCNext[control][target];

    // any single-qubit gate independent of the CNOT should be applied prior.
    for (std::size_t q = 0U; q < N; ++q) {
      if (q == control || q == target) {
        continue;
      }
      for (const auto& gate : SINGLE_QUBIT_GATES) {
        if (gate == qc::OpType::None) {
          continue;
        }
        disallowed = disallowed && !gSNext[gateToIndex(gate)][q];
      }
    }

    // no X gate may be placed on the target qubit since it would commute.
    disallowed = disallowed && !gSNext[gateToIndex(qc::OpType::X)][target];

    // no diagonal gate may be placed on the control qubit since it would
    // commute.
    disallowed = disallowed && !gSNext[gateToIndex(qc::OpType::Z)][control] &&
                 !gSNext[gateToIndex(qc::OpType::S)][control] &&
                 !gSNext[gateToIndex(qc::OpType::Sdg)][control];

    // no CNOT with the same control and a lower target qubit may be placed.
    for (std::size_t t = 0U; t < target; ++t) {
      if (t == control) {
        continue;
      }
      disallowed = disallowed && !gCNext[control][t];
    }

    // no CNOT with a lower control different from target may be placed.
    for (std::size_t c = 0U; c < control; ++c) {
      for (std::size_t t = 0U; t < N; ++t) {
        if (c == t || c == target) {
          continue;
        }
        disallowed = disallowed && !gCNext[c][t];
      }
    }

    lb->assertFormula(LogicTerm::implies(current, disallowed));
  }
}

} // namespace cs::encoding
