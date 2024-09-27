//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "ir/operations/OpType.hpp"
#include "logicblocks/Logic.hpp"
#include "logicblocks/LogicBlock.hpp"

#include <cstddef>
#include <memory>
#include <utility>

namespace cs::encoding {

class TableauEncoder {
public:
  TableauEncoder() = default;
  TableauEncoder(const std::size_t nQubits, const std::size_t tableauSize,
                 const std::size_t timestepLimit,
                 std::shared_ptr<logicbase::LogicBlock> logicBlock)
      : N(nQubits), S(tableauSize), T(timestepLimit),
        lb(std::move(logicBlock)) {}

  struct Variables {
    // variables for the X parts of the tableaus
    logicbase::LogicMatrix x;
    // variables for the Z parts of the tableaus
    logicbase::LogicMatrix z;
    // variables for the phase parts of the tableaus
    logicbase::LogicVector r;
    // variables for mapping of qubits
    logicbase::LogicMatrix p;

    // update rules for single-qubit gates
    [[nodiscard]] logicbase::LogicTerm
    singleQubitXChange(std::size_t pos, std::size_t qubit,
                       qc::OpType gate) const;
    [[nodiscard]] logicbase::LogicTerm
    singleQubitZChange(std::size_t pos, std::size_t qubit,
                       qc::OpType gate) const;
    [[nodiscard]] logicbase::LogicTerm
    singleQubitRChange(std::size_t pos, std::size_t qubit,
                       qc::OpType gate) const;

    // update rules for two-qubit gates
    [[nodiscard]] std::pair<logicbase::LogicTerm, logicbase::LogicTerm>
    twoQubitXChange(std::size_t pos, std::size_t ctrl, std::size_t trgt) const;
    [[nodiscard]] std::pair<logicbase::LogicTerm, logicbase::LogicTerm>
    twoQubitZChange(std::size_t pos, std::size_t ctrl, std::size_t trgt) const;
    [[nodiscard]] logicbase::LogicTerm
    twoQubitRChange(std::size_t pos, std::size_t ctrl, std::size_t trgt) const;
  };

  // variable creation
  void createTableauVariables();

  // fixing the tableau
  void assertTableau(const Tableau& tableau, std::size_t t) const;

  // extracting the tableau
  void extractTableauFromModel(Results& results, std::size_t t,
                               logicbase::Model& model) const;

  [[nodiscard]] auto* getVariables() { return &vars; }

  // get mapping variables and store them in results
  void extractMappingFromModel(Results& results, logicbase::Model& model) const;

  // assert constraints for mapping variables
  void assertMappingConstraints() const;

protected:
  // number of qubits N
  std::size_t N{}; // NOLINT (readability-identifier-naming)
  // number of rows in the tableau S
  std::size_t S{}; // NOLINT (readability-identifier-naming)
  // timestep limit T
  std::size_t T{}; // NOLINT (readability-identifier-naming)

  // the tableau variables
  Variables vars{};

  // the logic block to use
  std::shared_ptr<logicbase::LogicBlock> lb;
};
} // namespace cs::encoding
