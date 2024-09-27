//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"

#include "Logic.hpp"
#include "ir/operations/OpType.hpp"
#include "logicblocks/LogicTerm.hpp"
#include "sc/utils.hpp"

#include <cstddef>
#include <cstdint>
#include <plog/Log.h>
#include <plog/Severity.h>
#include <string>

namespace cs::encoding {

using namespace logicbase;

void encoding::MultiGateEncoder::assertConsistency() const {
  PLOG_DEBUG << "Asserting gate consistency";
  for (std::size_t t = 0U; t < T; ++t) {
    // asserting only a single gate is applied on each qubit.
    for (std::size_t q = 0U; q < N; ++q) {
      LogicVector gateVariables{};
      vars.collectSingleQubitGateVariables(t, q, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, true, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, false, gateVariables);

      IF_PLOG(plog::verbose) {
        PLOG_VERBOSE << "Gate variables at time " << t << " and qubit " << q;
        for (const auto& var : gateVariables) {
          PLOG_VERBOSE << var.getName();
        }
      }

      assertExactlyOne(gateVariables);
    }
  }
}

void encoding::MultiGateEncoder::assertGateConstraints() {
  PLOG_DEBUG << "Asserting gate constraints";
  xorHelpers = logicbase::LogicMatrix{T};
  for (std::size_t t = 0U; t < T; ++t) {
    PLOG_VERBOSE << "Asserting gate constraints at time " << t;
    rChanges = tvars->r[t + 1];
    splitXorR(tvars->r[t + 1], t);
    assertSingleQubitGateConstraints(t);
    assertTwoQubitGateConstraints(t);
    PLOG_VERBOSE << "Asserting r changes at time " << t;
    lb->assertFormula(tvars->r[t + 2] == xorHelpers[t].back());
  }
}

void encoding::MultiGateEncoder::assertSingleQubitGateConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    assertZConstraints(pos, q);
    assertXConstraints(pos, q);
    assertRConstraints(pos, q);
  }
}

void MultiGateEncoder::assertRConstraints(const std::size_t pos,
                                          const std::size_t qubit) {
  for (const auto gate : SINGLE_QUBIT_GATES) {
    const auto& change =
        LogicTerm::ite(vars.gS[pos][gateToIndex(gate)][qubit],
                       tvars->singleQubitRChange(pos + 1, qubit, gate),
                       LogicTerm(0, static_cast<std::uint16_t>(S)));
    splitXorR(change, pos);
  }
}

void encoding::MultiGateEncoder::assertTwoQubitGateConstraints(
    const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      // if no connection between ctrl and trgt then assert variable is false
      if (couplingMap.find(Edge{ctrl, trgt}) == couplingMap.end()) {
        PLOG_DEBUG << "Asserting no CNOT on " << ctrl << " and " << trgt;
        lb->assertFormula(!twoQubitGates[ctrl][trgt]);
        continue;
      }
      const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);
      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));

      PLOG_DEBUG << "Asserting CNOT on " << ctrl << " and " << trgt;
    }
  }
}

LogicTerm encoding::MultiGateEncoder::createTwoQubitGateConstraint(
    std::size_t pos, std::size_t ctrl, std::size_t trgt) {
  auto changes = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos + 1, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos + 1, ctrl, trgt);

  changes = changes && (tvars->x[pos + 2][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 2][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 2][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 2][trgt] == zTrgt);

  const auto& newRChanges = LogicTerm::ite(
      vars.gC[pos][ctrl][trgt], tvars->twoQubitRChange(pos + 1, ctrl, trgt),
      LogicTerm(0, static_cast<std::uint16_t>(S)));
  splitXorR(newRChanges, pos);
  return changes;
}

void MultiGateEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1];

  // once no gate is applied to the considered qubit, no single qubit gate can
  // be applied to it in the next time step.
  auto noSingleQubitGate = LogicTerm(true);
  for (const auto gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    noSingleQubitGate = noSingleQubitGate && !gSNext[gateToIndex(gate)][qubit];
  }
  lb->assertFormula(LogicTerm::implies(
      gSNow[gateToIndex(qc::OpType::None)][qubit], noSingleQubitGate));
}

void MultiGateEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow = vars.gS[pos];
  const auto& gCNext = vars.gC[pos + 1];

  // two identical CNOTs may not be applied in a row because they would cancel.
  lb->assertFormula(
      LogicTerm::implies(vars.gC[pos][ctrl][trgt], !gCNext[ctrl][trgt]));
  lb->assertFormula(
      LogicTerm::implies(vars.gC[pos][trgt][ctrl], !gCNext[trgt][ctrl]));

  // no gate on both qubits => no CNOT on them in the next time step.
  // hadamards on both qubits => no CNOT on them in the next time step (CNOT can
  // be conjugated with Hadamards) No Combination of Paulis on both Qubits
  // before a CNOT These gates can be just pushed through to the other side
  constexpr auto noneIndex = gateToIndex(qc::OpType::None);
  const auto noGate = gSNow[noneIndex][ctrl] && gSNow[noneIndex][trgt];
  const auto noFurtherCnot = !gCNext[ctrl][trgt] && !gCNext[trgt][ctrl];

  constexpr auto hIndex = gateToIndex(qc::OpType::H);
  constexpr auto xIndex = gateToIndex(qc::OpType::X);
  constexpr auto zIndex = gateToIndex(qc::OpType::Z);
  constexpr auto yIndex = gateToIndex(qc::OpType::Y);
  const auto hh = gSNow[hIndex][ctrl] && gSNow[hIndex][trgt];
  const auto gateBeforeCtrl =
      gSNow[zIndex][ctrl] || gSNow[xIndex][ctrl] || gSNow[yIndex][ctrl];
  const auto gateBeforeTarget =
      gSNow[zIndex][trgt] || gSNow[xIndex][trgt] || gSNow[yIndex][ctrl];
  lb->assertFormula(LogicTerm::implies(
      noGate || hh || (gateBeforeCtrl && gateBeforeTarget), noFurtherCnot));
}

void MultiGateEncoder::splitXorR(const logicbase::LogicTerm& changes,
                                 std::size_t pos) {
  auto& xorHelper = xorHelpers[pos];
  const std::string hName =
      "h_" + std::to_string(pos) + "_" + std::to_string(xorHelper.size());
  PLOG_DEBUG << "Creating helper variable for RChange XOR " << hName;
  const auto n = static_cast<std::uint16_t>(S);
  xorHelper.emplace_back(lb->makeVariable(hName, CType::BITVECTOR, n));
  if (xorHelper.size() == 1) {
    lb->assertFormula(xorHelper.back() == changes);
  } else {
    lb->assertFormula(xorHelper.back() ==
                      (xorHelper[xorHelpers[pos].size() - 2] ^ changes));
  }
}

} // namespace cs::encoding
