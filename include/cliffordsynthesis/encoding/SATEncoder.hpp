//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/Configuration.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/TargetMetric.hpp"
#include "cliffordsynthesis/encoding/GateEncoder.hpp"
#include "cliffordsynthesis/encoding/ObjectiveEncoder.hpp"
#include "cliffordsynthesis/encoding/TableauEncoder.hpp"
#include "logicblocks/Logic.hpp"
#include "logicblocks/LogicBlock.hpp"
#include "sc/utils.hpp"

#include <cstddef>
#include <memory>
#include <optional>

namespace cs::encoding {

class SATEncoder {
public:
  struct Configuration {
    // the initial tableau to encode at t=0
    Tableau* initialTableau{};

    // the target tableau to encode at t=T
    Tableau* targetTableau{};

    // the number of qubits to encode
    std::size_t nQubits{};

    // coupling Map of Qubits
    CouplingMap couplingMap;

    // the number of timesteps to encode
    std::size_t timestepLimit{};

    // minimal time steps of circuit
    std::size_t minimalTimesteps{};

    // the metric to optimize
    TargetMetric targetMetric = TargetMetric::Gates;

    // whether to use MaxSAT or Binary Search
    bool useMaxSAT = false;

    // whether to allow multiple gates per timestep
    bool useMultiGateEncoding = false;

    // whether to use symmetry breaking
    bool useSymmetryBreaking = false;

    // an optional limit on the total number of gates
    std::optional<std::size_t> gateLimit = std::nullopt;

    // an optional limit on the total number of two-qubit gates
    std::optional<std::size_t> twoQubitGateLimit = std::nullopt;

    SolverParameterMap solverParameters;
  };

  SATEncoder() = default;
  explicit SATEncoder(const Configuration& configuration)
      : config(configuration), N(configuration.nQubits),
        T(configuration.timestepLimit) {}
  virtual ~SATEncoder() = default;

  virtual Results run();

protected:
  void initializeSolver();
  void createFormulation();
  [[nodiscard]] logicbase::Result solve() const;
  void extractResultsFromModel(Results& res) const;
  void cleanup() const;

  std::shared_ptr<logicbase::LogicBlock> lb;
  std::shared_ptr<TableauEncoder> tableauEncoder;
  std::shared_ptr<GateEncoder> gateEncoder;
  std::shared_ptr<ObjectiveEncoder> objectiveEncoder;

  // all configuration options for the encoder
  Configuration config{};

  // number of qubits N
  std::size_t N{}; // NOLINT (readability-identifier-naming)
  // timestep limit T
  std::size_t T{}; // NOLINT (readability-identifier-naming)
};

} // namespace cs::encoding
