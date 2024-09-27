//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SATEncoder.hpp"

#include "Logic.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"
#include "cliffordsynthesis/encoding/ObjectiveEncoder.hpp"
#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"
#include "cliffordsynthesis/encoding/TableauEncoder.hpp"
#include "logicblocks/util_logicblock.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <plog/Log.h>
#include <stdexcept>
#include <string>
#include <variant>

namespace cs::encoding {

using namespace logicbase;

void SATEncoder::initializeSolver() {
  PLOG_DEBUG << "Initializing solver engine.";
  bool success = false;
  logicutil::Params params;
  for (const auto& [key, value] : config.solverParameters) {
    if (std::holds_alternative<bool>(value)) {
      params.addParam(key, std::get<bool>(value));
    } else if (std::holds_alternative<std::uint32_t>(value)) {
      params.addParam(key, std::get<std::uint32_t>(value));
    } else if (std::holds_alternative<std::string>(value)) {
      params.addParam(key, std::get<std::string>(value));
    } else if (std::holds_alternative<double>(value)) {
      params.addParam(key, std::get<double>(value));
    } else {
      const auto msg = "Unknown parameter type for " + key;
      PLOG_FATAL << msg;
      throw std::runtime_error(msg);
    }
  }

  if (config.useMaxSAT) {
    lb = logicutil::getZ3LogicOptimizer(success, true, params);
  } else {
    lb = logicutil::getZ3LogicBlock(success, true, params);
  }
  if (!success) {
    const auto* const msg = "Could not initialize solver engine.";
    PLOG_FATAL << msg;
    throw std::runtime_error(msg);
  }
}

void SATEncoder::createFormulation() {
  PLOG_INFO << "Creating formulation.";
  const auto start = std::chrono::high_resolution_clock::now();
  initializeSolver();

  const std::size_t s = config.targetTableau->hasDestabilizers() &&
                                config.initialTableau->hasDestabilizers()
                            ? 2U * N
                            : N;

  tableauEncoder = std::make_shared<TableauEncoder>(N, s, T + 2, lb);
  tableauEncoder->createTableauVariables();
  tableauEncoder->assertTableau(*config.initialTableau, 0U);
  tableauEncoder->assertTableau(*config.targetTableau, T + 2);
  tableauEncoder->assertMappingConstraints();

  if (config.useMultiGateEncoding) {
    gateEncoder = std::make_shared<MultiGateEncoder>(
        N, s, T, tableauEncoder->getVariables(), lb, config.couplingMap);
  } else {
    gateEncoder = std::make_shared<SingleGateEncoder>(
        N, s, T, tableauEncoder->getVariables(), lb, config.couplingMap);
  }
  gateEncoder->createSingleQubitGateVariables();
  gateEncoder->createTwoQubitGateVariables();
  gateEncoder->encodeGates();

  if (config.useSymmetryBreaking) {
    gateEncoder->encodeSymmetryBreakingConstraints();
  }

  objectiveEncoder =
      std::make_shared<ObjectiveEncoder>(N, T, gateEncoder->getVariables(), lb);

  if (config.gateLimit.has_value()) {
    objectiveEncoder->limitGateCount(*config.gateLimit, std::less_equal{});
  }

  if (config.twoQubitGateLimit.has_value()) {
    objectiveEncoder->limitGateCount(*config.twoQubitGateLimit,
                                     std::less_equal{}, false);
  }

  if (config.useMaxSAT) {
    objectiveEncoder->optimizeMetric(config.targetMetric);
  }

  const auto end = std::chrono::high_resolution_clock::now();
  const auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  PLOG_INFO << "Formulation created in " << duration << " ms.";
}

Result SATEncoder::solve() const {
  PLOG_INFO << "Solving the SAT instance.";

  const auto start = std::chrono::high_resolution_clock::now();
  const auto result = lb->solve();
  const auto end = std::chrono::high_resolution_clock::now();
  const auto runtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  PLOG_INFO << "Instance solved in " << runtime << " ms.";
  return result;
}

void SATEncoder::extractResultsFromModel(Results& res) const {
  auto* const model = lb->getModel();
  // result is tableau with mapping applied
  tableauEncoder->extractTableauFromModel(res, T + 1, *model);
  tableauEncoder->extractMappingFromModel(res, *model);
  gateEncoder->extractCircuitFromModel(res, *model);
}

void SATEncoder::cleanup() const {
  if (lb) {
    lb->reset();
  }
}
Results SATEncoder::run() {
  const auto start = std::chrono::high_resolution_clock::now();

  createFormulation();
  const auto solverResult = solve();

  const auto end = std::chrono::high_resolution_clock::now();
  const auto runtime = std::chrono::duration<double>(end - start);

  Results res{};
  res.setRuntime(runtime.count());
  res.setSolverResult(solverResult);

  if (solverResult == Result::SAT) {
    extractResultsFromModel(res);
  }

  cleanup();

  return res;
}

} // namespace cs::encoding
