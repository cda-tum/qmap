//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/encoding/TableauEncoder.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "logicblocks/Logic.hpp"
#include "logicblocks/LogicBlock.hpp"
#include "logicblocks/LogicTerm.hpp"
#include "sc/utils.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace cs::encoding {

class GateEncoder {
public:
  GateEncoder(const std::size_t nQubits, const std::size_t tableauSize,
              const std::size_t timestepLimit,
              TableauEncoder::Variables* tableauVars,
              std::shared_ptr<logicbase::LogicBlock> logicBlock, CouplingMap cm)
      : N(nQubits), S(tableauSize), T(timestepLimit),
        couplingMap(std::move(cm)), tvars(tableauVars),
        lb(std::move(logicBlock)) {}
  virtual ~GateEncoder() = default;

  struct Variables {
    // variables for the single-qubit gates
    logicbase::LogicMatrix3D gS;
    // variables for the two-qubit gates
    logicbase::LogicMatrix3D gC;

    void
    collectSingleQubitGateVariables(std::size_t pos, std::size_t qubit,
                                    logicbase::LogicVector& variables) const;
    void collectTwoQubitGateVariables(std::size_t pos, std::size_t qubit,
                                      bool target,
                                      logicbase::LogicVector& variables) const;
  };

  // variable creation
  void createSingleQubitGateVariables();
  void createTwoQubitGateVariables();

  // encode the relation between the tableaus and the gates
  virtual void encodeGates() {
    assertConsistency();
    assertGateConstraints();
  }

  virtual void encodeSymmetryBreakingConstraints();

  // extracting the circuit
  void extractCircuitFromModel(Results& res, logicbase::Model& model);

  [[nodiscard]] auto* getVariables() { return &vars; }

  static constexpr std::array<qc::OpType, 7> SINGLE_QUBIT_GATES = {
      qc::OpType::None, qc::OpType::X, qc::OpType::Y,  qc::OpType::Z,
      qc::OpType::H,    qc::OpType::S, qc::OpType::Sdg};

  [[nodiscard]] static constexpr std::size_t
  gateToIndex(const qc::OpType type) {
    for (std::size_t i = 0; i < SINGLE_QUBIT_GATES.size(); ++i) {
      if (SINGLE_QUBIT_GATES.at(i) == type) {
        return i;
      }
    }
    return 0;
  }

  template <qc::OpType Gate>
  [[nodiscard]] static constexpr bool containsGate() {
    for (const auto& g : // NOLINT(readability-use-anyofallof)
         SINGLE_QUBIT_GATES) {
      if (g == Gate) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] static constexpr bool containsX() {
    return containsGate<qc::OpType::X>();
  }
  [[nodiscard]] static constexpr bool containsY() {
    return containsGate<qc::OpType::Y>();
  }
  [[nodiscard]] static constexpr bool containsZ() {
    return containsGate<qc::OpType::Z>();
  }
  [[nodiscard]] static constexpr bool containsH() {
    return containsGate<qc::OpType::H>();
  }
  [[nodiscard]] static constexpr bool containsS() {
    return containsGate<qc::OpType::S>();
  }
  [[nodiscard]] static constexpr bool containsSdag() {
    return containsGate<qc::OpType::Sdg>();
  }

protected:
  // number of qubits N
  std::size_t N{}; // NOLINT (readability-identifier-naming)
  // number of rows in the tableau S
  std::size_t S{}; // NOLINT (readability-identifier-naming)
  // timestep limit T
  std::size_t T{}; // NOLINT (readability-identifier-naming)

  // coupling Map
  CouplingMap couplingMap;

  // the gate variables
  Variables vars{};

  // the tableau variables
  TableauEncoder::Variables* tvars{};

  // the logic block to use
  std::shared_ptr<logicbase::LogicBlock> lb;

  using TransformationFamily =
      std::pair<logicbase::LogicTerm, std::vector<qc::OpType>>;
  using GateToTransformation =
      std::function<logicbase::LogicTerm(std::size_t, std::size_t, qc::OpType)>;

  void assertExactlyOne(const logicbase::LogicVector& variables) const;

  virtual void assertConsistency() const = 0;

  virtual void assertGateConstraints() = 0;
  virtual void assertSingleQubitGateConstraints(std::size_t pos) = 0;
  virtual void assertTwoQubitGateConstraints(std::size_t pos) = 0;
  [[nodiscard]] static std::vector<TransformationFamily>
  collectGateTransformations(std::size_t pos, std::size_t qubit,
                             const GateToTransformation& gateToTransformation);
  void assertGatesImplyTransform(
      std::size_t pos, std::size_t qubit,
      const std::vector<TransformationFamily>& transformations);
  virtual void assertZConstraints(std::size_t pos, std::size_t qubit);
  virtual void assertXConstraints(std::size_t pos, std::size_t qubit);
  virtual void assertRConstraints(std::size_t pos, std::size_t qubit);
  [[nodiscard]] virtual logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) = 0;

  void extractSingleQubitGatesFromModel(std::size_t pos,
                                        logicbase::Model& model,
                                        qc::QuantumComputation& qc,
                                        std::size_t& nSingleQubitGates);
  void extractTwoQubitGatesFromModel(std::size_t pos, logicbase::Model& model,
                                     qc::QuantumComputation& qc,
                                     std::size_t& nTwoQubitGates);

  virtual void
  assertSingleQubitGateSymmetryBreakingConstraints(std::size_t pos);
  virtual void assertTwoQubitGateSymmetryBreakingConstraints(std::size_t pos);

  virtual void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                                     std::size_t qubit) = 0;
  virtual void assertSingleQubitGateCancellationConstraints(std::size_t pos,
                                                            std::size_t qubit);

  virtual void assertTwoQubitGateOrderConstraints(std::size_t pos,
                                                  std::size_t ctrl,
                                                  std::size_t trgt) = 0;
};
} // namespace cs::encoding
