//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/Tableau.hpp"
#include "ir/QuantumComputation.hpp"
#include "logicblocks/Logic.hpp"

#include <cstddef>
#include <limits>
#include <nlohmann/json_fwd.hpp>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cs {
class Results {
public:
  Results() = default;
  Results(qc::QuantumComputation& qc, const Tableau& tableau);

  virtual ~Results() = default;

  [[nodiscard]] std::size_t getGates() const {
    return getSingleQubitGates() + getTwoQubitGates();
  }
  [[nodiscard]] std::size_t getTwoQubitGates() const { return twoQubitGates; }
  [[nodiscard]] std::size_t getSingleQubitGates() const {
    return singleQubitGates;
  }
  [[nodiscard]] std::size_t getDepth() const { return depth; }
  [[nodiscard]] double getRuntime() const { return runtime; }
  [[nodiscard]] logicbase::Result getSolverResult() const {
    return solverResult;
  }
  [[nodiscard]] std::size_t getSolverCalls() const { return solverCalls; }

  [[nodiscard]] std::string getResultCircuit() const { return resultCircuit; }
  [[nodiscard]] std::string getResultTableau() const { return resultTableau; }

  [[nodiscard]] std::string getMapping() const {
    std::ostringstream oss;
    for (const auto& row : permutationVector) {
      for (const bool val : row) {
        oss << (val ? '1' : '0');
      }
      oss << '\n';
    }
    return oss.str();
  }
  [[nodiscard]] std::vector<std::vector<bool>> getMappingVector() const {
    return permutationVector;
  }

  void setSingleQubitGates(const std::size_t g) { singleQubitGates = g; }
  void setTwoQubitGates(const std::size_t g) { twoQubitGates = g; }
  void setDepth(const std::size_t d) { depth = d; }
  void setRuntime(const double t) { runtime = t; }
  void setSolverResult(const logicbase::Result r) { solverResult = r; }
  void setSolverCalls(const std::size_t c) { solverCalls = c; }

  void setResultCircuit(qc::QuantumComputation& qc);
  void setResultTableau(const Tableau& tableau);
  void setMapping(std::vector<std::vector<bool>> p) {
    std::ostringstream oss;
    for (const auto& row : permutationVector) {
      for (const bool val : row) {
        oss << (val ? '1' : '0');
      }
      oss << '\n';
    }
    permutationString = oss.str();
    permutationVector = std::move(p);
  }

  [[nodiscard]] bool sat() const {
    return getSolverResult() == logicbase::Result::SAT;
  }
  [[nodiscard]] bool unsat() const {
    return getSolverResult() == logicbase::Result::UNSAT;
  }

  [[nodiscard]] virtual nlohmann::basic_json<> json() const;

  friend std::ostream& operator<<(std::ostream& os, const Results& config);

protected:
  logicbase::Result solverResult = logicbase::Result::NDEF;
  std::size_t singleQubitGates = std::numeric_limits<std::size_t>::max();
  std::size_t twoQubitGates = std::numeric_limits<std::size_t>::max();
  std::size_t depth = std::numeric_limits<std::size_t>::max();
  double runtime = 0.0;
  std::size_t solverCalls = 0U;

  std::string permutationString;
  std::vector<std::vector<bool>> permutationVector;
  std::string resultTableau;
  std::string resultCircuit;
};

} // namespace cs
