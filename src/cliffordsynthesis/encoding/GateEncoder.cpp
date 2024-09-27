//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include "Definitions.hpp"
#include "Logic.hpp"
#include "LogicTerm.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "logicblocks/Encodings.hpp"
#include "logicblocks/Model.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <cstddef>
#include <plog/Log.h>
#include <string>
#include <vector>

namespace cs::encoding {

using namespace logicbase;

void GateEncoder::createSingleQubitGateVariables() {
  PLOG_DEBUG << "Creating single-qubit gate variables.";
  vars.gS.reserve(T);
  for (std::size_t t = 0U; t < T; ++t) {
    auto& timeStep = vars.gS.emplace_back();
    timeStep.reserve(SINGLE_QUBIT_GATES.size());
    for (const auto gate : SINGLE_QUBIT_GATES) {
      auto& g = timeStep.emplace_back();
      g.reserve(N);
      for (std::size_t q = 0U; q < N; ++q) {
        const std::string gName = "g_" + std::to_string(t) + "_" +
                                  toString(gate) + "_" + std::to_string(q);
        PLOG_VERBOSE << "Creating variable " << gName;
        g.emplace_back(lb->makeVariable(gName));
      }
    }
  }
}

void GateEncoder::createTwoQubitGateVariables() {
  PLOG_DEBUG << "Creating two-qubit gate variables.";
  vars.gC.reserve(T);
  for (std::size_t t = 0U; t < T; ++t) {
    auto& timeStep = vars.gC.emplace_back();
    timeStep.reserve(N);
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      auto& control = timeStep.emplace_back();
      control.reserve(N);
      for (std::size_t trgt = 0U; trgt < N; ++trgt) {
        const std::string gName = "g_" + std::to_string(t) + "_cx_" +
                                  std::to_string(ctrl) + "_" +
                                  std::to_string(trgt);
        PLOG_VERBOSE << "Creating variable " << gName;
        control.emplace_back(lb->makeVariable(gName));
      }
    }
  }
}

void GateEncoder::Variables::collectSingleQubitGateVariables(
    const std::size_t pos, const std::size_t qubit,
    LogicVector& variables) const {
  const auto& singleQubitGates = gS[pos];
  for (const auto& gate : singleQubitGates) {
    variables.emplace_back(gate[qubit]);
  }
}

void GateEncoder::Variables::collectTwoQubitGateVariables(
    const std::size_t pos, const std::size_t qubit, const bool target,
    LogicVector& variables) const {
  const auto& twoQubitGates = gC[pos];
  const auto n = twoQubitGates.size();
  for (std::size_t q = 0; q < n; ++q) {
    if (q == qubit) {
      continue;
    }
    if (target) {
      variables.emplace_back(twoQubitGates[q][qubit]);
    } else {
      variables.emplace_back(twoQubitGates[qubit][q]);
    }
  }
}

void GateEncoder::assertExactlyOne(const LogicVector& variables) const {
  const auto variableGrouping = encodings::groupVars(variables, 3U);
  lb->assertFormula(encodings::exactlyOneCmdr(variableGrouping,
                                              LogicTerm::noneTerm(), lb.get()));
}

std::vector<GateEncoder::TransformationFamily>
GateEncoder::collectGateTransformations(
    const std::size_t pos, const std::size_t qubit,
    const GateToTransformation& gateToTransformation) {
  std::vector<TransformationFamily> transformations;

  for (const auto& gate : SINGLE_QUBIT_GATES) {
    const auto& transformation = gateToTransformation(pos, qubit, gate);
    const auto& it = std::find_if(
        transformations.begin(), transformations.end(), [&](const auto& entry) {
          return entry.first.deepEquals(transformation);
        });
    if (it != transformations.end()) {
      it->second.emplace_back(gate);
    } else {
      transformations.emplace_back(transformation,
                                   std::vector<qc::OpType>{gate});
    }
  }
  return transformations;
}

void GateEncoder::assertGatesImplyTransform(
    const std::size_t pos, const std::size_t qubit,
    const std::vector<TransformationFamily>& transformations) {
  const auto& singleQubitGates = vars.gS[pos];
  for (const auto& [transformation, gates] : transformations) {
    auto gateOr = LogicTerm(false);
    for (const auto& gate : gates) {
      gateOr = gateOr || singleQubitGates[gateToIndex(gate)][qubit];
    }
    lb->assertFormula(LogicTerm::implies(gateOr, transformation));
  }
}

void GateEncoder::assertZConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToZTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitZChange(p1 + 1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToZTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->z[pos + 2][qubit] == transformation;
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::assertXConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToXTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitXChange(p1 + 1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToXTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->x[pos + 2][qubit] == transformation;
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::assertRConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToRTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitRChange(p1 + 1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToRTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->r[pos + 2] == (tvars->r[pos + 1] ^ transformation);
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates);
    extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates);
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(qc.getDepth());
  res.setResultCircuit(qc);
}

void GateEncoder::extractSingleQubitGatesFromModel(
    const std::size_t pos, Model& model, qc::QuantumComputation& qc,
    std::size_t& nSingleQubitGates) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      if (gate == qc::OpType::None) {
        continue;
      }
      if (model.getBoolValue(singleQubitGates[gateToIndex(gate)][q],
                             lb.get())) {
        qc.emplace_back<qc::StandardOperation>(q, gate);
        ++nSingleQubitGates;
        PLOG_DEBUG << toString(gate) << "(" << q << ")";
      }
    }
  }
}

void GateEncoder::extractTwoQubitGatesFromModel(const std::size_t pos,
                                                Model& model,
                                                qc::QuantumComputation& qc,
                                                size_t& nTwoQubitGates) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto control =
          qc::Control{static_cast<qc::Qubit>(ctrl), qc::Control::Type::Pos};
      if (model.getBoolValue(twoQubitGates[ctrl][trgt], lb.get())) {
        qc.cx(control, static_cast<qc::Qubit>(trgt));
        ++nTwoQubitGates;
        PLOG_DEBUG << "CX(" << ctrl << ", " << trgt << ")";
      }
    }
  }
}

void GateEncoder::encodeSymmetryBreakingConstraints() {
  PLOG_DEBUG << "Encoding symmetry breaking constraints.";
  for (std::size_t t = 0U; t < T; ++t) {
    assertSingleQubitGateSymmetryBreakingConstraints(t);
    assertTwoQubitGateSymmetryBreakingConstraints(t);
  }
}

void GateEncoder::assertSingleQubitGateCancellationConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert for the last timestep
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next timestep
  const auto& gSNow = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1U];

  // Any Pauli must not be followed by another Pauli since -iXYZ = I.
  std::vector<qc::OpType> paulis{};
  if constexpr (containsX()) {
    paulis.emplace_back(qc::OpType::X);
  }
  if constexpr (containsY()) {
    paulis.emplace_back(qc::OpType::Y);
  }
  if constexpr (containsZ()) {
    paulis.emplace_back(qc::OpType::Z);
  }
  constexpr bool containsPaulis = containsX() || containsY() || containsZ();
  if constexpr (containsPaulis) {
    auto gates = gSNow[gateToIndex(paulis[0])][qubit];
    auto disallowed = !gSNext[gateToIndex(paulis[0])][qubit];
    for (std::size_t i = 1U; i < paulis.size(); ++i) {
      gates = gates || gSNow[gateToIndex(paulis[i])][qubit];
      disallowed = disallowed && !gSNext[gateToIndex(paulis[i])][qubit];
    }

    if constexpr (containsH()) {
      // -(X|Y|Z)-H- ~= -H-(Z|Y|X)-
      constexpr auto gateIndex = gateToIndex(qc::OpType::H);
      disallowed = disallowed && !gSNext[gateIndex][qubit];
    }

    if constexpr (containsS() && containsSdag()) {
      const auto gateIndexS = gateToIndex(qc::OpType::S);
      const auto gateIndexSdg = gateToIndex(qc::OpType::Sdg);

      // -X-(S|Sd)- ~= -(Sd|S)-X-
      // -Y-(S|Sd)- ~= -(Sd|S)-Y-
      // -Z-(S|Sd)-  = -(S|Sd)-Z-
      disallowed = disallowed && !gSNext[gateIndexS][qubit] &&
                   !gSNext[gateIndexSdg][qubit];
    }

    lb->assertFormula(LogicTerm::implies(gates, disallowed));
  }

  // H is self-inverse
  if constexpr (containsH()) {
    constexpr auto gateIndex = gateToIndex(qc::OpType::H);
    lb->assertFormula(
        LogicTerm::implies(gSNow[gateIndex][qubit], !gSNext[gateIndex][qubit]));
  }

  if constexpr (containsS()) {
    constexpr auto gateIndexS = gateToIndex(qc::OpType::S);

    if constexpr (containsZ()) {
      constexpr auto gateIndexZ = gateToIndex(qc::OpType::Z);

      // -S-S- = -Z-
      auto gates = gSNow[gateIndexS][qubit];
      auto disallowed = !gSNext[gateIndexS][qubit];

      if constexpr (containsSdag()) {
        constexpr auto gateIndexSdag = gateToIndex(qc::OpType::Sdg);

        // -Sd-Sd- = -Z-
        // -Sd-S-  = -I-
        // -Sd-Z-  = -S-
        // -S-Sd-  = -I-
        // -S-Z-   = -Sd-
        gates = gates || gSNow[gateIndexSdag][qubit];
        disallowed = disallowed && !gSNext[gateIndexSdag][qubit] &&
                     !gSNext[gateIndexZ][qubit];
      }

      lb->assertFormula(LogicTerm::implies(gates, disallowed));
    } else {
      if constexpr (containsSdag()) {
        constexpr auto gateIndexSdag = gateToIndex(qc::OpType::Sdg);

        // -S-Sd- = -I-
        // -Sd-S- = -I-
        lb->assertFormula(LogicTerm::implies(gSNow[gateIndexS][qubit],
                                             !gSNext[gateIndexSdag][qubit]));
        lb->assertFormula(LogicTerm::implies(gSNow[gateIndexSdag][qubit],
                                             !gSNext[gateIndexS][qubit]));
      }
    }
  }
}

void GateEncoder::assertSingleQubitGateSymmetryBreakingConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    assertSingleQubitGateOrderConstraints(pos, q);
    assertSingleQubitGateCancellationConstraints(pos, q);
  }
}

void GateEncoder::assertTwoQubitGateSymmetryBreakingConstraints(
    const std::size_t pos) {
  for (std::size_t ctrl = 1U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < ctrl; ++trgt) {
      // avoid unnecessary constraints if CNOT cannot be applied anyways
      if (couplingMap.find(Edge{ctrl, trgt}) == couplingMap.end()) {
        continue;
      }
      assertTwoQubitGateOrderConstraints(pos, ctrl, trgt);
    }
  }
}

} // namespace cs::encoding
