//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "sc/heuristic/HeuristicMapper.hpp"

#include "Definitions.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "sc/Architecture.hpp"
#include "sc/DataLogger.hpp"
#include "sc/Mapper.hpp"
#include "sc/configuration/Configuration.hpp"
#include "sc/configuration/EarlyTermination.hpp"
#include "sc/configuration/Heuristic.hpp"
#include "sc/configuration/InitialLayout.hpp"
#include "sc/configuration/Layering.hpp"
#include "sc/configuration/LookaheadHeuristic.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <utility>
#include <vector>

void HeuristicMapper::map(const Configuration& configuration) {
  if (configuration.dataLoggingEnabled()) {
    dataLogger = std::make_unique<DataLogger>(configuration.dataLoggingPath,
                                              *architecture, qc);
  }

  tightHeur = isTight(configuration.heuristic);
  fidelityAwareHeur = isFidelityAware(configuration.heuristic);

  results = MappingResults{};
  results.config = configuration;
  const auto& config = results.config;
  checkParameters();
  const auto start = std::chrono::steady_clock::now();
  initResults();

  // perform pre-mapping optimizations
  preMappingOptimizations(config);

  createLayers();
  if (config.verbose) {
    std::clog << "Teleportation qubits: " << config.teleportationQubits << "\n";
    printLayering(std::clog);
  }

  createInitialMapping();
  if (config.verbose) {
    printLocations(std::clog);
    printQubits(std::clog);
  }

  for (std::size_t i = 0; i < config.iterativeBidirectionalRoutingPasses; ++i) {
    if (config.verbose) {
      std::clog << "\nIterative bidirectional routing (forward pass " << i
                << "):\n";
    }
    pseudoRouteCircuit(false);
    if (config.verbose) {
      std::clog << "\nIterative bidirectional routing (backward pass " << i
                << "):\n";
    }
    pseudoRouteCircuit(true);

    if (config.verbose) {
      std::clog << "\nMain routing:\n";
    }
  }

  routeCircuit();

  postMappingOptimizations(config);
  countGates(qcMapped, results.output);
  finalizeMappedCircuit();

  const auto end = std::chrono::steady_clock::now();
  const std::chrono::duration<double> diff = end - start;
  results.time = diff.count();
  results.timeout = false;

  if (config.dataLoggingEnabled()) {
    dataLogger->logOutputCircuit(qcMapped);
    dataLogger->logMappingResult(results);
    dataLogger->close();
  }
}

void HeuristicMapper::staticInitialMapping() {
  for (const auto& gate : layers.at(0U)) {
    if (gate.singleQubit()) {
      continue;
    }

    for (const auto& [q0, q1] : architecture->getCouplingMap()) {
      if (qubits.at(q0) == DEFAULT_POSITION &&
          qubits.at(q1) == DEFAULT_POSITION) {
        qubits.at(q0) = gate.control;
        qubits.at(q1) = static_cast<std::int16_t>(gate.target);
        locations.at(static_cast<std::uint16_t>(gate.control)) =
            static_cast<std::int16_t>(q0);
        locations.at(gate.target) = static_cast<std::int16_t>(q1);
        qcMapped.initialLayout.at(q0) = static_cast<qc::Qubit>(gate.control);
        qcMapped.initialLayout.at(q1) = static_cast<qc::Qubit>(gate.target);
        qcMapped.outputPermutation.at(q0) =
            static_cast<qc::Qubit>(gate.control);
        qcMapped.outputPermutation.at(q1) = static_cast<qc::Qubit>(gate.target);
        break;
      }
    }
  }

  // assign remaining logical qubits
  for (qc::Qubit i = 0U; i < architecture->getNqubits(); ++i) {
    if (qc.initialLayout.count(i) > 0 && locations.at(i) == DEFAULT_POSITION) {
      for (qc::Qubit j = 0U; j < architecture->getNqubits(); ++j) {
        if (qubits.at(j) == DEFAULT_POSITION) {
          locations.at(i) = static_cast<std::int16_t>(j);
          qubits.at(j) = static_cast<std::int16_t>(i);
          qcMapped.initialLayout.at(j) = i;
          qcMapped.outputPermutation.at(j) = i;
          break;
        }
      }
    }
  }
}

void HeuristicMapper::checkParameters() {
  const auto& config = results.config;
  if (config.layering == Layering::OddGates ||
      config.layering == Layering::QubitTriangle) {
    throw QMAPException("Layering strategy " + toString(config.layering) +
                        " not suitable for heuristic mapper!");
  }
  if (fidelityAwareHeur && !architecture->isFidelityAvailable()) {
    throw QMAPException("Fidelity aware heuristic chosen, but no or "
                        "insufficient calibration data available for this "
                        "architecture!");
  }
  if (fidelityAwareHeur && !isFidelityAware(config.lookaheadHeuristic) &&
      config.lookaheadHeuristic != LookaheadHeuristic::None) {
    throw QMAPException("Fidelity-aware heuristics may only be used with "
                        "fidelity-aware lookahead heuristics (or no "
                        "lookahead)!");
  }
  if (fidelityAwareHeur && config.teleportationQubits > 0) {
    throw QMAPException("Teleportation is not yet supported for heuristic "
                        "mapper using fidelity-aware mapping!");
  }
}

void HeuristicMapper::createInitialMapping() {
  auto& config = results.config;

  if (layers.empty()) {
    return;
  }

  if (config.teleportationQubits > 0) {
    std::mt19937_64 mt;
    if (config.teleportationSeed == 0) {
      std::array<std::mt19937_64::result_type, std::mt19937_64::state_size>
          randomData{};
      std::random_device rd;
      std::generate(std::begin(randomData), std::end(randomData),
                    [&rd]() { return rd(); });
      std::seed_seq seeds(std::begin(randomData), std::end(randomData));
      mt.seed(seeds);
    } else {
      mt.seed(config.teleportationSeed);
    }

    std::uniform_int_distribution<> dis(
        0,
        static_cast<std::int32_t>(architecture->getCouplingMap().size() - 1));

    for (std::size_t i = 0; i < config.teleportationQubits; i += 2) {
      Edge e{};
      do { // NOLINT(cppcoreguidelines-avoid-do-while)
        auto it = std::begin(architecture->getCouplingMap());
        std::advance(it, dis(mt));
        e = *it;
      } while (qubits.at(e.first) != -1 || qubits.at(e.second) != -1);
      const auto teleportationQubit = qc.getNqubits() + i;
      locations.at(teleportationQubit) = static_cast<std::int16_t>(e.first);
      locations.at(teleportationQubit + 1) =
          static_cast<std::int16_t>(e.second);
      qubits.at(e.first) = static_cast<std::int16_t>(teleportationQubit);
      qubits.at(e.second) = static_cast<std::int16_t>(teleportationQubit + 1);
    }

    if (config.teleportationFake) {
      config.teleportationQubits = 0;
    }
  }

  switch (config.initialLayout) {
  case InitialLayout::Identity:
    for (qc::Qubit i = 0; i < architecture->getNqubits(); ++i) {
      if (qc.initialLayout.count(i) > 0) {
        locations.at(i) = static_cast<std::int16_t>(i);
        qubits.at(i) = static_cast<std::int16_t>(i);
      }
    }
    break;
  case InitialLayout::Static:
    staticInitialMapping();
    break;
  case InitialLayout::Dynamic:
    // nothing to be done here
    break;

    // TODO: Design strategy that maps most used qubit to most connected qubits
    // on architecture
  }
}

namespace {
// search for current position of target value in map and afterward exchange
// it with the value at new position
void findAndSWAP(const qc::Qubit targetValue, const qc::Qubit newPosition,
                 qc::Permutation& map) {
  for (const auto& q : map) {
    if (q.second == targetValue) {
      std::swap(map.at(newPosition), map.at(q.first));
      break;
    }
  }
}
} // namespace

void HeuristicMapper::mapUnmappedGates(std::size_t layer) {
  if (fidelityAwareHeur) {
    for (std::size_t q = 0; q < singleQubitMultiplicities.at(layer).size();
         ++q) {
      if (singleQubitMultiplicities.at(layer).at(q) == 0) {
        continue;
      }
      if (locations.at(q) == DEFAULT_POSITION) {
        // TODO: consider fidelity
        // map to first free physical qubit
        for (std::uint16_t physQbit = 0; physQbit < architecture->getNqubits();
             ++physQbit) {
          if (qubits.at(physQbit) == -1) {
            locations.at(q) = static_cast<std::int16_t>(physQbit);
            qubits.at(physQbit) = static_cast<std::int16_t>(q);
            break;
          }
        }
      }
    }
  }

  for (const auto& [logEdge, _] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = logEdge;

    const auto q1Location = locations.at(q1);
    const auto q2Location = locations.at(q2);

    if (q1Location == DEFAULT_POSITION && q2Location == DEFAULT_POSITION) {
      std::set<Edge> possibleEdges{};
      // gather all edges in the architecture for which both qubits are unmapped
      for (const auto& edge : architecture->getCouplingMap()) {
        if (qubits.at(edge.first) == DEFAULT_POSITION &&
            qubits.at(edge.second) == DEFAULT_POSITION) {
          possibleEdges.emplace(edge);
        }
      }
      std::pair<std::uint16_t, std::uint16_t> chosenEdge;

      if (possibleEdges.empty()) {
        // map to 2 qubits with minimal distance
        double bestScore = std::numeric_limits<int>::max();

        for (std::uint16_t i = 0; i < architecture->getNqubits(); i++) {
          for (std::uint16_t j = i + 1; j < architecture->getNqubits(); j++) {
            if (qubits.at(i) == DEFAULT_POSITION &&
                qubits.at(j) == DEFAULT_POSITION) {
              const double dist = architecture->distance(i, j);
              if (dist < bestScore) {
                bestScore = dist;
                chosenEdge = std::make_pair(i, j);
              }
            }
          }
        }
      } else {
        chosenEdge = *possibleEdges.begin();
      }
      // TODO: Consider fidelity here if available. The best available edge
      // should be chosen
      locations.at(q1) = static_cast<std::int16_t>(chosenEdge.first);
      locations.at(q2) = static_cast<std::int16_t>(chosenEdge.second);
      qubits.at(chosenEdge.first) = static_cast<std::int16_t>(q1);
      qubits.at(chosenEdge.second) = static_cast<std::int16_t>(q2);
      findAndSWAP(q1, chosenEdge.first, qcMapped.initialLayout);
      findAndSWAP(q2, chosenEdge.second, qcMapped.initialLayout);
      findAndSWAP(q1, chosenEdge.first, qcMapped.outputPermutation);
      findAndSWAP(q2, chosenEdge.second, qcMapped.outputPermutation);
    } else if (q1Location == DEFAULT_POSITION) {
      mapToMinDistance(q2, q1);
    } else if (q2Location == DEFAULT_POSITION) {
      mapToMinDistance(q1, q2);
    }
  }
}

void HeuristicMapper::mapToMinDistance(const std::uint16_t source,
                                       const std::uint16_t target) {
  auto min = std::numeric_limits<double>::max();
  std::optional<std::uint16_t> pos = std::nullopt;
  for (std::uint16_t i = 0; i < architecture->getNqubits(); ++i) {
    if (qubits.at(i) == DEFAULT_POSITION) {
      // TODO: Consider fidelity here if available
      const auto distance = architecture->distance(
          static_cast<std::uint16_t>(locations.at(source)), i);
      if (distance < min) {
        min = distance;
        pos = i;
      }
    }
  }
  assert(pos.has_value());
  qubits.at(*pos) = static_cast<std::int16_t>(target);
  locations.at(target) = static_cast<std::int16_t>(*pos);
  findAndSWAP(target, *pos, qcMapped.initialLayout);
  findAndSWAP(target, *pos, qcMapped.outputPermutation);
}

void HeuristicMapper::pseudoRouteCircuit(bool reverse) {
  // save original global data for restoring it later
  const auto originalResults = results;
  const auto originalLayers = layers;
  const auto originalSingleQubitMultiplicities = singleQubitMultiplicities;
  const auto originalTwoQubitMultiplicities = twoQubitMultiplicities;
  const auto originalActiveQubits = activeQubits;
  const auto originalActiveQubits1QGates = activeQubits1QGates;
  const auto originalActiveQubits2QGates = activeQubits2QGates;

  auto& config = results.config;
  config.dataLoggingPath = ""; // disable data logging for pseudo routing
  config.debug = false;

  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto layerIndex = (reverse ? layers.size() - i - 1 : i);
    const Node result = aStarMap(layerIndex, reverse);

    qubits = result.qubits;
    locations = result.locations;

    if (config.verbose) {
      printLocations(std::clog);
      printQubits(std::clog);
    }
  }

  // restore original global data
  results = originalResults;
  layers = originalLayers;
  singleQubitMultiplicities = originalSingleQubitMultiplicities;
  twoQubitMultiplicities = originalTwoQubitMultiplicities;
  activeQubits = originalActiveQubits;
  activeQubits1QGates = originalActiveQubits1QGates;
  activeQubits2QGates = originalActiveQubits2QGates;
}

void HeuristicMapper::routeCircuit() {
  const auto& config = results.config;

  std::size_t gateidx = 0;
  std::vector<std::size_t> gatesToAdjust{};
  results.output.gates = 0U;
  for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
    const Node result = aStarMap(layerIndex, false);

    qubits = result.qubits;
    locations = result.locations;

    if (config.verbose) {
      printLocations(std::clog);
      printQubits(std::clog);
    }

    if (layerIndex != 0 && config.addBarriersBetweenLayers) {
      qcMapped.barrier();
      gateidx++;
    }

    // initial layer needs no swaps
    if (layerIndex != 0 || config.swapOnFirstLayer) {
      for (const auto& swap : result.swaps) {
        if (swap.op == qc::SWAP) {
          if (config.verbose) {
            std::clog << "SWAP: " << swap.first << " <-> " << swap.second
                      << "\n";
          }
          // check if SWAP is placed on a valid edge
          assert(
              architecture->isEdgeConnected({swap.first, swap.second}, false));
          qcMapped.swap(swap.first, swap.second);
          results.output.swaps++;
        } else if (swap.op == qc::Teleportation) {
          if (config.verbose) {
            std::clog << "TELE: " << swap.first << " <-> " << swap.second
                      << "\n";
          }
          qcMapped.emplace_back<qc::StandardOperation>(
              qc::Targets{static_cast<qc::Qubit>(swap.first),
                          static_cast<qc::Qubit>(swap.second),
                          static_cast<qc::Qubit>(swap.middleAncilla)},
              qc::Teleportation);
          results.output.teleportations++;
        }
        gateidx++;
      }
    }

    // add gates of the layer to circuit
    for (const auto& gate : layers.at(layerIndex)) {
      auto* op = dynamic_cast<qc::StandardOperation*>(gate.op);
      if (op == nullptr) {
        throw QMAPException(
            "Cast to StandardOperation not possible during mapping. Check that "
            "circuit contains only StandardOperations");
      }

      if (gate.singleQubit()) {
        if (locations.at(gate.target) == DEFAULT_POSITION) {
          qcMapped.emplace_back<qc::StandardOperation>(
              gate.target, op->getType(), op->getParameter());
          gatesToAdjust.push_back(gateidx);
          gateidx++;
        } else {
          qcMapped.emplace_back<qc::StandardOperation>(
              locations.at(gate.target), op->getType(), op->getParameter());
          gateidx++;
        }
      } else {
        const Edge cnot = {
            locations.at(static_cast<std::uint16_t>(gate.control)),
            locations.at(gate.target)};
        if (!architecture->isEdgeConnected(cnot)) {
          const Edge reversed = {cnot.second, cnot.first};
          // check if CNOT is placed on a valid edge
          assert(architecture->isEdgeConnected(reversed));
          qcMapped.h(reversed.first);
          qcMapped.h(reversed.second);
          qcMapped.cx(qc::Control{static_cast<qc::Qubit>(reversed.first)},
                      reversed.second);
          qcMapped.h(reversed.second);
          qcMapped.h(reversed.first);

          results.output.directionReverse++;
          gateidx += 5;
        } else {
          qcMapped.cx(qc::Control{static_cast<qc::Qubit>(cnot.first)},
                      cnot.second);
          gateidx++;
        }
      }
    }
  }

  if (config.debug && results.heuristicBenchmark.expandedNodes > 0) {
    auto& benchmark = results.heuristicBenchmark;
    benchmark.secondsPerNode /= static_cast<double>(benchmark.expandedNodes);
    benchmark.averageBranchingFactor =
        static_cast<double>(benchmark.generatedNodes - layers.size()) /
        static_cast<double>(benchmark.expandedNodes);
    for (const auto& layer : results.layerHeuristicBenchmark) {
      benchmark.effectiveBranchingFactor +=
          layer.effectiveBranchingFactor *
          (static_cast<double>(layer.expandedNodes) /
           static_cast<double>(benchmark.expandedNodes));
    }
  }

  // infer output permutation from qubit locations
  qcMapped.outputPermutation.clear();
  for (std::size_t i = 0U; i < architecture->getNqubits(); ++i) {
    if (const auto lq = qubits.at(i); lq != -1) {
      const auto logicalQubit = static_cast<qc::Qubit>(lq);
      // check whether this is a qubit from the original circuit
      if (logicalQubit < qc.getNqubits()) {
        qcMapped.outputPermutation[static_cast<qc::Qubit>(i)] =
            static_cast<qc::Qubit>(qubits.at(i));
      } else {
        qcMapped.setLogicalQubitGarbage(logicalQubit);
      }
    }
  }

  // fix single qubit gates
  if (!gatesToAdjust.empty()) {
    gateidx--; // index of last operation
    for (auto it = qcMapped.rbegin(); it != qcMapped.rend(); ++it, --gateidx) {
      auto* op = dynamic_cast<qc::StandardOperation*>(it->get());
      if (op == nullptr) {
        throw QMAPException(
            "Cast to StandardOperation not possible during mapping. Check that "
            "circuit contains only StandardOperations");
      }
      if (op->getType() == qc::SWAP) {
        const auto q0 = qubits.at(op->getTargets().at(0));
        const auto q1 = qubits.at(op->getTargets().at(1));
        qubits.at(op->getTargets().at(0)) = q1;
        qubits.at(op->getTargets().at(1)) = q0;

        if (q0 != DEFAULT_POSITION) {
          locations.at(static_cast<std::size_t>(q0)) =
              static_cast<std::int16_t>(op->getTargets().at(1));
        }
        if (q1 != DEFAULT_POSITION) {
          locations.at(static_cast<std::size_t>(q1)) =
              static_cast<std::int16_t>(op->getTargets().at(0));
        }
      }
      if (!gatesToAdjust.empty() && gatesToAdjust.back() == gateidx) {
        gatesToAdjust.pop_back();
        const auto target = op->getTargets().at(0);
        const auto targetLocation = locations.at(target);

        if (targetLocation == -1) {
          // qubit only occurs in single qubit gates, can be mapped to an
          // arbitrary free qubit
          std::uint16_t loc = 0;
          while (qubits.at(loc) != DEFAULT_POSITION) {
            ++loc;
          }
          locations.at(target) = static_cast<std::int16_t>(loc);
          qubits.at(loc) = static_cast<std::int16_t>(target);
          op->setTargets({static_cast<qc::Qubit>(loc)});
          qcMapped.initialLayout.at(loc) = target;
          qcMapped.outputPermutation[static_cast<qc::Qubit>(loc)] = target;
        } else {
          op->setTargets({static_cast<qc::Qubit>(targetLocation)});
        }
      }
    }
  }

  // mark every qubit that is not mapped to a logical qubit as garbage
  std::size_t count = 0U;
  for (std::size_t i = 0U; i < architecture->getNqubits(); ++i) {
    if (const auto lq = qubits.at(i); lq == -1) {
      qcMapped.setLogicalQubitGarbage(
          static_cast<qc::Qubit>(qc.getNqubits() + count));
      ++count;
    }
  }
}

HeuristicMapper::Node HeuristicMapper::aStarMap(size_t layer, bool reverse) {
  const auto& config = results.config;
  nextNodeId = 0;

  const SingleQubitMultiplicity& singleQubitMultiplicity =
      singleQubitMultiplicities.at(layer);
  const TwoQubitMultiplicity& twoQubitMultiplicity =
      twoQubitMultiplicities.at(layer);
  Node node(architecture->getNqubits(), nextNodeId++);
  Node bestDoneNode(architecture->getNqubits(), 0);
  bool validMapping = false;

  mapUnmappedGates(layer);

  node.locations = locations;
  node.qubits = qubits;
  recalculateFixedCost(layer, node);
  updateHeuristicCost(layer, node);
  updateLookaheadPenalty(layer, node);

  if (config.dataLoggingEnabled()) {
    dataLogger->logSearchNode(layer, node.id, node.parent,
                              node.costFixed + node.costFixedReversals,
                              node.costHeur, node.lookaheadPenalty, node.qubits,
                              node.validMapping, node.swaps, node.depth);
  }
  nodes.push(node);

  const auto start = std::chrono::steady_clock::now();
  std::size_t expandedNodes = 0;
  std::size_t expandedNodesAfterFirstSolution = 0;
  std::size_t expandedNodesAfterOptimalSolution = 0;
  std::size_t solutionNodes = 0;
  std::size_t solutionNodesAfterOptimalSolution = 0;
  bool earlyTermination = false;

  const bool splittable =
      config.automaticLayerSplits ? isLayerSplittable(layer) : false;

  while (!nodes.empty() &&
         (!validMapping ||
          nodes.top().getTotalCost() < bestDoneNode.getTotalFixedCost())) {
    if (splittable && expandedNodes >= config.automaticLayerSplitsNodeLimit) {
      if (config.dataLoggingEnabled()) {
        qc::CompoundOperation compOp{};
        for (const auto& gate : layers.at(layer)) {
          compOp.emplace_back(gate.op->clone());
        }

        dataLogger->logFinalizeLayer(layer, compOp, singleQubitMultiplicity,
                                     twoQubitMultiplicity, qubits, 0, 0, 0, 0,
                                     {}, {}, 0);
        dataLogger->splitLayer();
      }
      splitLayer(layer, *architecture);
      if (config.verbose) {
        std::clog << "Split layer\n";
      }
      // recursively restart search with newly split layer
      // (step to the end of the circuit, if reverse mapping is active, since
      // the split layer is inserted in this direction, otherwise 1 layer would
      // be skipped)
      return aStarMap(reverse ? layer + 1 : layer, reverse);
    }
    Node current = nodes.top();
    if (current.validMapping) {
      ++solutionNodes;
      if (!validMapping ||
          current.getTotalFixedCost() < bestDoneNode.getTotalFixedCost()) {
        bestDoneNode = current;
        expandedNodesAfterOptimalSolution = 0;
        solutionNodesAfterOptimalSolution = 0;
      } else {
        ++solutionNodesAfterOptimalSolution;
      }
      validMapping = true;
      if (tightHeur) {
        break;
      }
    }
    nodes.pop();
    expandNode(current, layer);
    ++expandedNodes;
    if (validMapping) {
      ++expandedNodesAfterFirstSolution;
      ++expandedNodesAfterOptimalSolution;

      if (config.earlyTermination != EarlyTermination::None) {
        std::size_t n = 0;
        switch (config.earlyTermination) {
        case EarlyTermination::ExpandedNodes:
          n = expandedNodes;
          break;
        case EarlyTermination::ExpandedNodesAfterFirstSolution:
          n = expandedNodesAfterFirstSolution;
          break;
        case EarlyTermination::ExpandedNodesAfterCurrentOptimalSolution:
          n = expandedNodesAfterOptimalSolution;
          break;
        case EarlyTermination::SolutionNodes:
          n = solutionNodes;
          break;
        case EarlyTermination::SolutionNodesAfterCurrentOptimalSolution:
          n = solutionNodesAfterOptimalSolution;
          break;
        default:
          break;
        }
        if (n >= config.earlyTerminationLimit) {
          earlyTermination = true;
          break;
        }
      }
    }
  }

  if (!validMapping) {
    throw QMAPException("No viable mapping found.");
  }

  Node result = bestDoneNode;
  if (config.debug) {
    const auto end = std::chrono::steady_clock::now();
    results.layerHeuristicBenchmark.emplace_back();
    auto layerResultsIt = results.layerHeuristicBenchmark.rbegin();
    layerResultsIt->expandedNodes = expandedNodes;
    results.heuristicBenchmark.expandedNodes += expandedNodes;

    layerResultsIt->solutionDepth = result.depth;
    layerResultsIt->expandedNodesAfterFirstSolution =
        expandedNodesAfterFirstSolution;
    layerResultsIt->expandedNodesAfterOptimalSolution =
        expandedNodesAfterOptimalSolution;
    layerResultsIt->solutionNodes = solutionNodes;
    layerResultsIt->solutionNodesAfterOptimalSolution =
        solutionNodesAfterOptimalSolution;
    layerResultsIt->earlyTermination = earlyTermination;

    const std::chrono::duration<double> diff = end - start;
    results.heuristicBenchmark.secondsPerNode += diff.count();

    layerResultsIt->generatedNodes = nextNodeId;
    results.heuristicBenchmark.generatedNodes += layerResultsIt->generatedNodes;

    if (layerResultsIt->expandedNodes > 0) {
      layerResultsIt->secondsPerNode =
          diff.count() / static_cast<double>(layerResultsIt->expandedNodes);
      layerResultsIt->averageBranchingFactor =
          static_cast<double>(layerResultsIt->generatedNodes - 1) /
          static_cast<double>(layerResultsIt->expandedNodes);
    }

    layerResultsIt->effectiveBranchingFactor = computeEffectiveBranchingRate(
        layerResultsIt->expandedNodes + 1, result.depth);
  }

  if (config.dataLoggingEnabled()) {
    qc::CompoundOperation compOp{};
    for (const auto& gate : layers.at(layer)) {
      compOp.emplace_back(gate.op->clone());
    }

    dataLogger->logFinalizeLayer(
        layer, compOp, singleQubitMultiplicities.at(layer),
        twoQubitMultiplicities.at(layer), qubits, result.id, result.costFixed,
        result.costHeur, result.lookaheadPenalty, result.qubits, result.swaps,
        result.depth);
  }

  // clear nodes
  while (!nodes.empty()) {
    nodes.pop();
  }

  return result;
}

void HeuristicMapper::expandNode(Node& node, std::size_t layer) {
  const auto& consideredQubits = getConsideredQubits(layer);
  std::vector<std::vector<bool>> usedSwaps;
  usedSwaps.reserve(architecture->getNqubits());
  for (int p = 0; p < architecture->getNqubits(); ++p) {
    usedSwaps.emplace_back(architecture->getNqubits());
  }

  // set up new teleportation qubits
  std::set<Edge> perms = architecture->getCouplingMap();
  architecture->getCurrentTeleportations().clear();
  architecture->getTeleportationQubits().clear();
  for (std::size_t i = 0; i < results.config.teleportationQubits; i += 2) {
    architecture->getTeleportationQubits().emplace_back(
        node.locations.at(qc.getNqubits() + i),
        node.locations.at(qc.getNqubits() + i + 1));
    Edge e;
    for (auto const& g : architecture->getCouplingMap()) {
      if (g.first == node.locations.at(qc.getNqubits() + i) &&
          g.second != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first = g.second;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i) &&
          g.first != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first = g.first;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.first == node.locations.at(qc.getNqubits() + i + 1) &&
          g.second != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.second;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i + 1) &&
          g.first != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.first;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
    }
  }

  for (const auto& q : consideredQubits) {
    for (const auto& edge : perms) {
      if (edge.first == node.locations.at(q) ||
          edge.second == node.locations.at(q)) {
        const auto q1 = node.qubits.at(edge.first);
        const auto q2 = node.qubits.at(edge.second);
        if (q2 == -1 || q1 == -1) {
          expandNodeAddOneSwap(edge, node, layer);
        } else if (!usedSwaps.at(static_cast<std::size_t>(q1))
                        .at(static_cast<std::size_t>(q2))) {
          usedSwaps.at(static_cast<std::size_t>(q1))
              .at(static_cast<std::size_t>(q2)) = true;
          usedSwaps.at(static_cast<std::size_t>(q2))
              .at(static_cast<std::size_t>(q1)) = true;
          expandNodeAddOneSwap(edge, node, layer);
        }
      }
    }
  }
}

void HeuristicMapper::expandNodeAddOneSwap(const Edge& swap, Node& node,
                                           const std::size_t layer) {
  Node newNode =
      Node(nextNodeId++, node.id, node.qubits, node.locations, node.swaps,
           node.validMappedTwoQubitGates, node.costFixed,
           node.costFixedReversals, node.depth + 1, node.sharedSwaps);

  if (architecture->isEdgeConnected(swap, false)) {
    applySWAP(swap, layer, newNode);
  } else {
    applyTeleportation(swap, layer, newNode);
  }

  nodes.push(newNode);
  if (results.config.dataLoggingEnabled()) {
    dataLogger->logSearchNode(layer, newNode.id, newNode.parent,
                              newNode.costFixed + newNode.costFixedReversals,
                              newNode.costHeur, newNode.lookaheadPenalty,
                              newNode.qubits, newNode.validMapping,
                              newNode.swaps, newNode.depth);
  }
}

void HeuristicMapper::recalculateFixedCost(std::size_t layer, Node& node) {
  node.validMappedTwoQubitGates.clear();
  for (const auto& [edge, mult] : twoQubitMultiplicities.at(layer)) {
    const auto [q1, q2] = edge;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    if (architecture->isEdgeConnected({physQ1, physQ2}, false)) {
      // validly mapped
      node.validMappedTwoQubitGates.emplace(q1, q2);
    }
  }

  if (fidelityAwareHeur) {
    recalculateFixedCostFidelity(layer, node);
  } else {
    recalculateFixedCostNonFidelity(node);
  }
  recalculateFixedCostReversals(layer, node);
}

void HeuristicMapper::recalculateFixedCostReversals(std::size_t layer,
                                                    Node& node) {
  node.costFixedReversals = 0.;
  if (architecture->bidirectional() || fidelityAwareHeur ||
      node.validMappedTwoQubitGates.size() !=
          twoQubitMultiplicities.at(layer).size()) {
    // costFixedReversals should only be non-zero in goal nodes for
    // non-fidelity-aware heuristics and if there are unidirectional
    // edges in the architecture
    return;
  }

  // only consider reversal costs as fixed in goal nodes
  for (const auto& [edge, mult] : twoQubitMultiplicities.at(layer)) {
    const auto [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = mult;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    if (!architecture->isEdgeConnected({physQ1, physQ2})) {
      node.costFixedReversals += forwardMult * COST_DIRECTION_REVERSE;
    } else if (!architecture->isEdgeConnected({physQ2, physQ1})) {
      node.costFixedReversals += reverseMult * COST_DIRECTION_REVERSE;
    }
  }
}

void HeuristicMapper::recalculateFixedCostNonFidelity(Node& node) {
  node.costFixed = 0;

  // swap costs
  for (auto& swap : node.swaps) {
    if (swap.op == qc::SWAP) {
      // branch clone intended for performance reasons (checking edge-wise for
      // bidirectionality is not O(1))
      // NOLINTBEGIN(bugprone-branch-clone)
      if (architecture->bidirectional()) {
        node.costFixed += COST_BIDIRECTIONAL_SWAP;
      } else if (architecture->unidirectional() ||
                 !architecture->isEdgeBidirectional(
                     {swap.first, swap.second})) {
        node.costFixed += COST_UNIDIRECTIONAL_SWAP;
      } else {
        node.costFixed += COST_BIDIRECTIONAL_SWAP;
      }
      // NOLINTEND(bugprone-branch-clone)
    } else if (swap.op == qc::Teleportation) {
      node.costFixed += COST_TELEPORTATION;
    }
  }
}

void HeuristicMapper::recalculateFixedCostFidelity(std::size_t layer,
                                                   Node& node) {
  const auto& singleQubitGateMultiplicity = singleQubitMultiplicities.at(layer);
  const auto& twoQubitGateMultiplicity = twoQubitMultiplicities.at(layer);

  node.costFixed = 0;
  // adding costs of single qubit gates
  for (std::uint16_t i = 0U; i < architecture->getNqubits(); ++i) {
    if (singleQubitGateMultiplicity.at(i) == 0) {
      continue;
    }
    node.costFixed += singleQubitGateMultiplicity.at(i) *
                      architecture->getSingleQubitFidelityCost(
                          static_cast<std::uint16_t>(node.locations.at(i)));
  }
  // adding cost of the swap gates
  for (auto& swap : node.swaps) {
    if (swap.op == qc::SWAP) {
      node.costFixed +=
          architecture->getSwapFidelityCost(swap.first, swap.second);
    } else if (swap.op == qc::Teleportation) {
      throw QMAPException("Teleportation currently not supported for "
                          "noise-aware mapping");
    }
  }
  // adding cost of two qubit gates that are already mapped next to each other
  for (const auto& [edge, mult] : twoQubitGateMultiplicity) {
    if (node.validMappedTwoQubitGates.find(edge) ==
        node.validMappedTwoQubitGates.end()) {
      // 2-qubit-gates not yet validly mapped are handled in the heuristic
      continue;
    }
    const auto [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = mult;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    node.costFixed +=
        (forwardMult * architecture->getTwoQubitFidelityCost(physQ1, physQ2) +
         reverseMult * architecture->getTwoQubitFidelityCost(physQ2, physQ1));
  }
}

void HeuristicMapper::applySWAP(const Edge& swap, std::size_t layer,
                                Node& node) {
  assert(architecture->isEdgeConnected(swap, false));
  const auto& singleQubitGateMultiplicity = singleQubitMultiplicities.at(layer);

  const auto q1 = node.qubits.at(swap.first);
  const auto q2 = node.qubits.at(swap.second);

  updateSharedSwaps(swap, layer, node);

  node.qubits.at(swap.first) = q2;
  node.qubits.at(swap.second) = q1;

  if (q1 != -1) {
    node.locations.at(static_cast<std::size_t>(q1)) =
        static_cast<std::int16_t>(swap.second);
  }
  if (q2 != -1) {
    node.locations.at(static_cast<std::size_t>(q2)) =
        static_cast<std::int16_t>(swap.first);
  }

  node.swaps.emplace_back(swap.first, swap.second, qc::SWAP);

  // check if swap created or destroyed any valid mappings of qubit pairs
  for (const auto& [edge, mult] : twoQubitMultiplicities.at(layer)) {
    const auto [q3, q4] = edge;
    if (q3 == q1 || q3 == q2 || q4 == q1 || q4 == q2) {
      const auto physQ3 = static_cast<std::uint16_t>(node.locations.at(q3));
      const auto physQ4 = static_cast<std::uint16_t>(node.locations.at(q4));
      if (architecture->isEdgeConnected({physQ3, physQ4}, false)) {
        // validly mapped now
        if (fidelityAwareHeur && node.validMappedTwoQubitGates.find(edge) ==
                                     node.validMappedTwoQubitGates.end()) {
          // not mapped validly before
          // add cost of newly validly mapped gates
          node.costFixed +=
              mult.first *
                  architecture->getTwoQubitFidelityCost(physQ3, physQ4) +
              mult.second *
                  architecture->getTwoQubitFidelityCost(physQ4, physQ3);
        }
        node.validMappedTwoQubitGates.emplace(edge);
      } else {
        // not mapped validly now
        if (fidelityAwareHeur && node.validMappedTwoQubitGates.find(edge) !=
                                     node.validMappedTwoQubitGates.end()) {
          // mapped validly before
          // remove cost of now no longer validly mapped gates
          auto prevPhysQ3 = physQ3;
          if (prevPhysQ3 == swap.first) {
            prevPhysQ3 = swap.second;
          } else if (prevPhysQ3 == swap.second) {
            prevPhysQ3 = swap.first;
          }
          auto prevPhysQ4 = physQ4;
          if (prevPhysQ4 == swap.first) {
            prevPhysQ4 = swap.second;
          } else if (prevPhysQ4 == swap.second) {
            prevPhysQ4 = swap.first;
          }

          node.costFixed -=
              mult.first * architecture->getTwoQubitFidelityCost(prevPhysQ3,
                                                                 prevPhysQ4) +
              mult.second *
                  architecture->getTwoQubitFidelityCost(prevPhysQ4, prevPhysQ3);
        }
        node.validMappedTwoQubitGates.erase(edge);
      }
    }
  }

  if (fidelityAwareHeur) {
    std::uint16_t q1Mult = 0;
    std::uint16_t q2Mult = 0;
    if (q1 != -1) {
      q1Mult = singleQubitGateMultiplicity.at(static_cast<std::size_t>(q1));
    }
    if (q2 != -1) {
      q2Mult = singleQubitGateMultiplicity.at(static_cast<std::size_t>(q2));
    }
    // accounting for fidelity difference of single qubit gates (two qubit
    // gates are handled in the heuristic)
    node.costFixed +=
        ((q2Mult - q1Mult) *
             architecture->getSingleQubitFidelityCost(swap.first) +
         (q1Mult - q2Mult) *
             architecture->getSingleQubitFidelityCost(swap.second));
    // adding cost of the swap gate itself
    node.costFixed +=
        architecture->getSwapFidelityCost(swap.first, swap.second);
  } else {
    // branch clone intended for performance reasons (checking edge-wise for
    // bidirectionality is not O(1))
    // NOLINTBEGIN(bugprone-branch-clone)
    if (architecture->bidirectional()) {
      node.costFixed += COST_BIDIRECTIONAL_SWAP;
    } else if (architecture->unidirectional() ||
               !architecture->isEdgeBidirectional({swap.first, swap.second})) {
      node.costFixed += COST_UNIDIRECTIONAL_SWAP;
    } else {
      node.costFixed += COST_BIDIRECTIONAL_SWAP;
    }
    // NOLINTEND(bugprone-branch-clone)
  }

  recalculateFixedCostReversals(layer, node);
  updateHeuristicCost(layer, node);
  if (results.config.lookaheadHeuristic != LookaheadHeuristic::None) {
    updateLookaheadPenalty(layer, node);
  }
}

void HeuristicMapper::applyTeleportation(const Edge& swap, std::size_t layer,
                                         Node& node) {
  const auto q1 = node.qubits.at(swap.first);
  const auto q2 = node.qubits.at(swap.second);

  updateSharedSwaps(swap, layer, node);

  node.qubits.at(swap.first) = q2;
  node.qubits.at(swap.second) = q1;

  if (q1 != -1) {
    node.locations.at(static_cast<std::size_t>(q1)) =
        static_cast<std::int16_t>(swap.second);
  }
  if (q2 != -1) {
    node.locations.at(static_cast<std::size_t>(q2)) =
        static_cast<std::int16_t>(swap.first);
  }

  std::uint16_t middleAnc = std::numeric_limits<decltype(middleAnc)>::max();
  for (const auto& qpair : architecture->getTeleportationQubits()) {
    if (swap.first == qpair.first || swap.second == qpair.first) {
      middleAnc = static_cast<std::uint16_t>(qpair.second);
    } else if (swap.first == qpair.second || swap.second == qpair.second) {
      middleAnc = static_cast<std::uint16_t>(qpair.first);
    }
  }

  // Teleportation between wrong qubits
  assert(middleAnc != std::numeric_limits<decltype(middleAnc)>::max());

  std::uint16_t source = std::numeric_limits<decltype(source)>::max();
  std::uint16_t target = std::numeric_limits<decltype(target)>::max();
  if (architecture->isEdgeConnected({swap.first, middleAnc}, false)) {
    source = swap.first;
    target = swap.second;
  } else {
    source = swap.second;
    target = swap.first;
  }

  if (source == middleAnc || target == middleAnc) {
    std::clog << "FAIL: TELE " << source << " -(" << middleAnc << ")-> "
              << target << "\n";
    throw QMAPException("Overlap between source/target and middle "
                        "ancillary in teleportation.");
  }

  node.swaps.emplace_back(source, target, middleAnc, qc::Teleportation);

  node.costFixed += COST_TELEPORTATION;

  // check if swap created or destroyed any valid mappings of qubit pairs
  for (const auto& [edge, mult] : twoQubitMultiplicities.at(layer)) {
    const auto [q3, q4] = edge;
    if (q3 == q1 || q3 == q2 || q4 == q1 || q4 == q2) {
      const auto physQ3 = static_cast<std::uint16_t>(node.locations.at(q3));
      const auto physQ4 = static_cast<std::uint16_t>(node.locations.at(q4));
      if (architecture->isEdgeConnected({physQ3, physQ4}, false)) {
        // validly mapped now
        node.validMappedTwoQubitGates.emplace(edge);
      } else {
        // not mapped validly now
        node.validMappedTwoQubitGates.erase(edge);
      }
    }
  }

  recalculateFixedCostReversals(layer, node);
  updateHeuristicCost(layer, node);
  if (results.config.lookaheadHeuristic != LookaheadHeuristic::None) {
    updateLookaheadPenalty(layer, node);
  }
}

void HeuristicMapper::updateSharedSwaps(const Edge& swap, std::size_t layer,
                                        Node& node) {
  const auto& consideredQubits = getConsideredQubits(layer);
  const auto& twoQubitGateMultiplicity = twoQubitMultiplicities.at(layer);

  const auto q1 = node.qubits.at(swap.first);
  const auto q2 = node.qubits.at(swap.second);
  if (q1 == -1 || q2 == -1 ||
      consideredQubits.find(static_cast<std::uint16_t>(q1)) ==
          consideredQubits.end() ||
      consideredQubits.find(static_cast<std::uint16_t>(q2)) ==
          consideredQubits.end()) {
    // the given swap can only be a shared swap if both qubits are active in
    // the current layer
    return;
  }

  // TODO: handle single qubit gates for fidelity aware heuristic, if
  //        `Node::sharedSwaps` is ever used in a fidelity aware heuristic
  Edge logEdge1 = {q1, q1};
  Edge logEdge2 = {q2, q2};
  for (const auto& [edge, multiplicity] : twoQubitGateMultiplicity) {
    if (edge.first == q1) {
      logEdge1.second = edge.second;
    } else if (edge.second == q1) {
      logEdge1.second = edge.first;
    }
    if (edge.first == q2) {
      logEdge2.second = edge.second;
    } else if (edge.second == q2) {
      logEdge2.second = edge.first;
    }
  }
  if ( // if both swapped qubits are acted on by a 2q gate
      logEdge1.second != q1 && logEdge2.second != q2 &&
      // if it is not the same 2q gate acting on both qubits
      logEdge1.second != q2) {
    auto physQ3 =
        static_cast<std::uint16_t>(node.locations.at(logEdge1.second));
    auto physQ4 =
        static_cast<std::uint16_t>(node.locations.at(logEdge2.second));

    double logEdge1DistanceBefore = 0.;
    double logEdge1DistanceNew = 0.;
    double logEdge2DistanceBefore = 0.;
    double logEdge2DistanceNew = 0.;
    if (fidelityAwareHeur) {
      logEdge1DistanceBefore =
          std::min(architecture->fidelityDistance(swap.first, physQ3),
                   architecture->fidelityDistance(physQ3, swap.first));
      logEdge1DistanceNew =
          std::min(architecture->fidelityDistance(swap.second, physQ3),
                   architecture->fidelityDistance(physQ3, swap.second));
      logEdge2DistanceBefore =
          std::min(architecture->fidelityDistance(swap.second, physQ4),
                   architecture->fidelityDistance(physQ4, swap.second));
      logEdge2DistanceNew =
          std::min(architecture->fidelityDistance(swap.first, physQ4),
                   architecture->fidelityDistance(physQ4, swap.first));
    } else {
      logEdge1DistanceBefore =
          std::min(architecture->distance(swap.first, physQ3, false),
                   architecture->distance(physQ3, swap.first, false));
      logEdge1DistanceNew =
          std::min(architecture->distance(swap.second, physQ3, false),
                   architecture->distance(physQ3, swap.second, false));
      logEdge2DistanceBefore =
          std::min(architecture->distance(swap.second, physQ4, false),
                   architecture->distance(physQ4, swap.second, false));
      logEdge2DistanceNew =
          std::min(architecture->distance(swap.first, physQ4, false),
                   architecture->distance(physQ4, swap.first, false));
    }
    if (logEdge1DistanceNew < logEdge1DistanceBefore &&
        logEdge2DistanceNew < logEdge2DistanceBefore) {
      ++node.sharedSwaps;
    }
  }
}

void HeuristicMapper::updateHeuristicCost(std::size_t layer, Node& node) {
  // the mapping is valid, only if all qubit pairs are mapped next to each other
  node.validMapping = (node.validMappedTwoQubitGates.size() ==
                       twoQubitMultiplicities.at(layer).size());

  switch (results.config.heuristic) {
  case Heuristic::GateCountMaxDistance:
    node.costHeur = heuristicGateCountMaxDistance(layer, node);
    break;
  case Heuristic::GateCountSumDistance:
    node.costHeur = heuristicGateCountSumDistance(layer, node);
    break;
  case Heuristic::GateCountSumDistanceMinusSharedSwaps:
    node.costHeur = heuristicGateCountSumDistanceMinusSharedSwaps(layer, node);
    break;
  case Heuristic::GateCountMaxDistanceOrSumDistanceMinusSharedSwaps:
    node.costHeur =
        heuristicGateCountMaxDistanceOrSumDistanceMinusSharedSwaps(layer, node);
    break;
  case Heuristic::FidelityBestLocation:
    node.costHeur = heuristicFidelityBestLocation(layer, node);
    break;
  default:
    throw QMAPException("Unknown heuristic.");
  }
}

double HeuristicMapper::heuristicGateCountMaxDistance(std::size_t layer,
                                                      Node& node) {
  if (node.validMapping) {
    return 0.;
  }
  double costHeur = 0.;

  for (const auto& [edge, multiplicity] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = multiplicity;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    if (!architecture->bidirectional() &&
        node.validMappedTwoQubitGates.find(edge) !=
            node.validMappedTwoQubitGates.end()) {
      // validly mapped 2-qubit-gates
      if (!architecture->isEdgeConnected({physQ1, physQ2})) {
        costHeur =
            std::max(costHeur,
                     static_cast<double>(forwardMult * COST_DIRECTION_REVERSE));
      } else if (!architecture->isEdgeConnected({physQ2, physQ1})) {
        costHeur =
            std::max(costHeur,
                     static_cast<double>(reverseMult * COST_DIRECTION_REVERSE));
      }
    } else {
      // not validly mapped 2-qubit-gates
      if (forwardMult > 0) {
        costHeur = std::max(costHeur, architecture->distance(physQ1, physQ2));
      }
      if (reverseMult > 0) {
        costHeur = std::max(costHeur, architecture->distance(physQ2, physQ1));
      }
    }
  }

  return costHeur;
}

double HeuristicMapper::heuristicGateCountSumDistance(std::size_t layer,
                                                      Node& node) {
  if (node.validMapping) {
    return 0.;
  }
  double costHeur = 0.;

  for (const auto& [edge, multiplicity] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = multiplicity;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    if (!architecture->bidirectional() &&
        node.validMappedTwoQubitGates.find(edge) !=
            node.validMappedTwoQubitGates.end()) {
      // validly mapped 2-qubit-gates
      if (!architecture->isEdgeConnected({physQ1, physQ2})) {
        costHeur += forwardMult * COST_DIRECTION_REVERSE;
      } else if (!architecture->isEdgeConnected({physQ2, physQ1})) {
        costHeur += reverseMult * COST_DIRECTION_REVERSE;
      }
    } else {
      // not validly mapped 2-qubit-gates
      double swapCost = 0.;

      if (forwardMult == 0) {
        // forwardMult == 0 && reverseMult > 0
        swapCost = architecture->distance(physQ2, physQ1);
      } else if (reverseMult == 0) {
        // forwardMult > 0 && reverseMult == 0
        swapCost = architecture->distance(physQ1, physQ2);
      } else {
        // forwardMult > 0 && reverseMult > 0
        swapCost = std::max(architecture->distance(physQ1, physQ2),
                            architecture->distance(physQ2, physQ1));
      }
      costHeur += swapCost;
    }
  }

  return costHeur;
}

double HeuristicMapper::heuristicGateCountSumDistanceMinusSharedSwaps(
    std::size_t layer, Node& node) {
  if (node.validMapping) {
    return 0.;
  }
  const auto& twoQubitGateMultiplicity = twoQubitMultiplicities.at(layer);
  double costHeur = 0.;
  double costReversals = 0.;
  std::vector<std::size_t> nSwaps{};
  nSwaps.reserve(twoQubitGateMultiplicity.size());

  for (const auto& [edge, multiplicity] : twoQubitGateMultiplicity) {
    const auto& [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = multiplicity;
    const auto physQ1 = static_cast<std::uint16_t>(node.locations.at(q1));
    const auto physQ2 = static_cast<std::uint16_t>(node.locations.at(q2));

    if (architecture->unidirectional()) {
      // only for purely unidirectional architectures is it certain that at
      // least one of the two directions has to be reversed
      costReversals +=
          std::min(forwardMult, reverseMult) * COST_DIRECTION_REVERSE;
    }

    if (node.validMappedTwoQubitGates.find(edge) !=
        node.validMappedTwoQubitGates.end()) {
      // validly mapped 2-qubit-gates
      continue;
    }

    double swapCost = 0.;
    if (forwardMult == 0) {
      // forwardMult == 0 && reverseMult > 0
      swapCost = architecture->distance(physQ2, physQ1, false);
    } else if (reverseMult == 0) {
      // forwardMult > 0 && reverseMult == 0
      swapCost = architecture->distance(physQ1, physQ2, false);
    } else {
      // forwardMult > 0 && reverseMult > 0
      swapCost = std::min(architecture->distance(physQ1, physQ2, false),
                          architecture->distance(physQ2, physQ1, false));
    }
    costHeur += swapCost;

    // infer maximum number of swaps in this distance
    if (architecture->unidirectional()) {
      nSwaps.emplace_back(
          static_cast<std::size_t>(swapCost / COST_UNIDIRECTIONAL_SWAP));
    } else {
      nSwaps.emplace_back(
          static_cast<std::size_t>(swapCost / COST_BIDIRECTIONAL_SWAP));
    }
  }

  // sort number of swaps in descending order
  std::sort(nSwaps.begin(), nSwaps.end(), std::greater<>());

  // infer maximum number of shared swaps
  std::size_t maxSharedSwaps = 0;
  for (std::size_t i = 0; i < nSwaps.size() - 1; ++i) {
    std::size_t maxSharedSwapsEdge =
        0; // maximum number of shared swaps for this edge
    for (std::size_t j = i + 1;
         j < nSwaps.size() && maxSharedSwapsEdge < nSwaps[i]; ++j) {
      if (nSwaps[j] > 0) {
        ++maxSharedSwapsEdge;
        --nSwaps[j];
      }
    }
    maxSharedSwaps += maxSharedSwapsEdge;
  }
  if (node.sharedSwaps < maxSharedSwaps) {
    maxSharedSwaps -= node.sharedSwaps;
  } else {
    maxSharedSwaps = 0;
  }

  double sharedSwapCostReduction = 0;
  if (architecture->bidirectional()) {
    sharedSwapCostReduction =
        static_cast<double>(maxSharedSwaps * COST_BIDIRECTIONAL_SWAP);
  } else {
    sharedSwapCostReduction =
        static_cast<double>(maxSharedSwaps * COST_UNIDIRECTIONAL_SWAP);
  }

  return std::max(0., costHeur - sharedSwapCostReduction) + costReversals;
}

double
HeuristicMapper::heuristicGateCountMaxDistanceOrSumDistanceMinusSharedSwaps(
    std::size_t layer, Node& node) {
  return std::max(heuristicGateCountMaxDistance(layer, node),
                  heuristicGateCountSumDistanceMinusSharedSwaps(layer, node));
}

double HeuristicMapper::heuristicFidelityBestLocation(std::size_t layer,
                                                      Node& node) {
  const auto& consideredQubits = getConsideredQubits(layer);
  const auto& singleQubitGateMultiplicity = singleQubitMultiplicities.at(layer);
  const auto& twoQubitGateMultiplicity = twoQubitMultiplicities.at(layer);

  double costHeur = 0.;

  // single qubit gate savings potential by moving them to different physical
  // qubits with higher fidelity
  double savingsPotential = 0.;
  for (std::uint16_t logQbit = 0U; logQbit < architecture->getNqubits();
       ++logQbit) {
    if (singleQubitGateMultiplicity.at(logQbit) == 0) {
      continue;
    }
    double qbitSavings = 0;
    const double currFidelity = architecture->getSingleQubitFidelityCost(
        static_cast<std::uint16_t>(node.locations.at(logQbit)));
    for (std::uint16_t physQbit = 0U; physQbit < architecture->getNqubits();
         ++physQbit) {
      if (architecture->getSingleQubitFidelityCost(physQbit) >= currFidelity) {
        continue;
      }
      const double curSavings =
          singleQubitGateMultiplicity.at(logQbit) *
              (currFidelity -
               architecture->getSingleQubitFidelityCost(physQbit)) -
          architecture->fidelityDistance(
              static_cast<std::uint16_t>(node.locations.at(logQbit)), physQbit,
              consideredQubits.size() - 1);
      qbitSavings = std::max(qbitSavings, curSavings);
    }
    savingsPotential += qbitSavings;
  }

  // iterating over all virtual qubit pairs, that share a gate on the
  // current layer
  for (const auto& [edge, mult] : twoQubitGateMultiplicity) {
    const auto [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = mult;

    const bool edgeDone = (node.validMappedTwoQubitGates.find(edge) !=
                               node.validMappedTwoQubitGates.end() ||
                           node.validMappedTwoQubitGates.find({q1, q2}) !=
                               node.validMappedTwoQubitGates.end());

    // find the optimal edge, to which to remap the given virtual qubit
    // pair and take the cost of moving it there via swaps plus the
    // fidelity cost  of executing all their shared gates on that edge
    // as the qubit pairs cost
    double swapCost = std::numeric_limits<double>::max();
    for (const auto& [q3, q4] : architecture->getCouplingMap()) {
      swapCost = std::min(
          swapCost,
          forwardMult * architecture->getTwoQubitFidelityCost(q3, q4) +
              reverseMult * architecture->getTwoQubitFidelityCost(q4, q3) +
              architecture->fidelityDistance(
                  static_cast<std::uint16_t>(node.locations.at(q1)), q3,
                  consideredQubits.size() - 1) +
              architecture->fidelityDistance(
                  static_cast<std::uint16_t>(node.locations.at(q2)), q4,
                  consideredQubits.size() - 1));
      swapCost = std::min(
          swapCost,
          forwardMult * architecture->getTwoQubitFidelityCost(q4, q3) +
              reverseMult * architecture->getTwoQubitFidelityCost(q3, q4) +
              architecture->fidelityDistance(
                  static_cast<std::uint16_t>(node.locations.at(q2)), q3,
                  consideredQubits.size() - 1) +
              architecture->fidelityDistance(
                  static_cast<std::uint16_t>(node.locations.at(q1)), q4,
                  consideredQubits.size() - 1));
    }

    if (edgeDone) {
      const double currEdgeCost =
          (forwardMult *
               architecture->getTwoQubitFidelityCost(
                   static_cast<std::uint16_t>(node.locations.at(q1)),
                   static_cast<std::uint16_t>(node.locations.at(q2))) +
           reverseMult *
               architecture->getTwoQubitFidelityCost(
                   static_cast<std::uint16_t>(node.locations.at(q2)),
                   static_cast<std::uint16_t>(node.locations.at(q1))));
      savingsPotential += (currEdgeCost - swapCost);
    } else {
      costHeur += swapCost;
    }
  }

  return costHeur - savingsPotential;
}

void HeuristicMapper::updateLookaheadPenalty(const std::size_t layer,
                                             HeuristicMapper::Node& node) {
  const auto& config = results.config;
  node.lookaheadPenalty = 0.;
  auto nextLayer = getNextLayer(layer);
  double factor = config.firstLookaheadFactor;

  for (std::size_t i = 0; i < config.nrLookaheads; ++i) {
    if (nextLayer == std::numeric_limits<std::size_t>::max()) {
      break;
    }

    double penalty = 0.;
    switch (config.lookaheadHeuristic) {
    case LookaheadHeuristic::GateCountMaxDistance:
      penalty = lookaheadGateCountMaxDistance(nextLayer, node);
      break;
    case LookaheadHeuristic::GateCountSumDistance:
      penalty = lookaheadGateCountSumDistance(nextLayer, node);
      break;
    default:
      break;
    }

    node.lookaheadPenalty += factor * penalty;
    factor *= config.lookaheadFactor;
    nextLayer = getNextLayer(nextLayer); // TODO: consider single qubits here
                                         // for better fidelity lookahead
  }
}

double
HeuristicMapper::lookaheadGateCountMaxDistance(const std::size_t layer,
                                               HeuristicMapper::Node& node) {
  double penalty = 0.;

  for (const auto& [edge, multiplicity] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = multiplicity;

    const auto loc1 = node.locations.at(q1);
    const auto loc2 = node.locations.at(q2);
    if (loc1 == DEFAULT_POSITION && loc2 == DEFAULT_POSITION) {
      // no penalty
    } else if (loc1 == DEFAULT_POSITION) {
      auto min = std::numeric_limits<double>::max();
      for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
        if (node.qubits.at(j) == DEFAULT_POSITION) {
          // TODO: Consider fidelity here if available
          if (forwardMult > 0) {
            min = std::min(min, architecture->distance(
                                    j, static_cast<std::uint16_t>(loc2)));
          }
          if (reverseMult > 0) {
            min = std::min(min, architecture->distance(
                                    static_cast<std::uint16_t>(loc2), j));
          }
        }
      }
      penalty = std::max(penalty, min);
    } else if (loc2 == DEFAULT_POSITION) {
      auto min = std::numeric_limits<double>::max();
      for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
        if (node.qubits.at(j) == DEFAULT_POSITION) {
          // TODO: Consider fidelity here if available
          if (forwardMult > 0) {
            min = std::min(min, architecture->distance(
                                    static_cast<std::uint16_t>(loc1), j));
          }
          if (reverseMult > 0) {
            min = std::min(min, architecture->distance(
                                    j, static_cast<std::uint16_t>(loc1)));
          }
        }
      }
      penalty = std::max(penalty, min);
    } else {
      double cost = std::numeric_limits<double>::max();
      if (forwardMult > 0) {
        cost = std::min(
            cost, architecture->distance(static_cast<std::uint16_t>(loc1),
                                         static_cast<std::uint16_t>(loc2)));
      }
      if (reverseMult > 0) {
        cost = std::min(
            cost, architecture->distance(static_cast<std::uint16_t>(loc2),
                                         static_cast<std::uint16_t>(loc1)));
      }
      penalty = std::max(penalty, cost);
    }
  }

  return penalty;
}

double
HeuristicMapper::lookaheadGateCountSumDistance(const std::size_t layer,
                                               HeuristicMapper::Node& node) {
  double penalty = 0.;

  for (const auto& [edge, multiplicity] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = edge;
    const auto [forwardMult, reverseMult] = multiplicity;

    const auto loc1 = node.locations.at(q1);
    const auto loc2 = node.locations.at(q2);
    if (loc1 == DEFAULT_POSITION && loc2 == DEFAULT_POSITION) {
      // no penalty
    } else if (loc1 == DEFAULT_POSITION) {
      auto min = std::numeric_limits<double>::max();
      for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
        if (node.qubits.at(j) == DEFAULT_POSITION) {
          // TODO: Consider fidelity here if available
          if (forwardMult > 0) {
            min = std::min(min, architecture->distance(
                                    j, static_cast<std::uint16_t>(loc2)));
          }
          if (reverseMult > 0) {
            min = std::min(min, architecture->distance(
                                    static_cast<std::uint16_t>(loc2), j));
          }
        }
      }
      penalty += min;
    } else if (loc2 == DEFAULT_POSITION) {
      auto min = std::numeric_limits<double>::max();
      for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
        if (node.qubits.at(j) == DEFAULT_POSITION) {
          // TODO: Consider fidelity here if available
          if (forwardMult > 0) {
            min = std::min(min, architecture->distance(
                                    static_cast<std::uint16_t>(loc1), j));
          }
          if (reverseMult > 0) {
            min = std::min(min, architecture->distance(
                                    j, static_cast<std::uint16_t>(loc1)));
          }
        }
      }
      penalty += min;
    } else {
      double cost = std::numeric_limits<double>::max();
      if (forwardMult > 0) {
        cost = std::min(
            cost, architecture->distance(static_cast<std::uint16_t>(loc1),
                                         static_cast<std::uint16_t>(loc2)));
      }
      if (reverseMult > 0) {
        cost = std::min(
            cost, architecture->distance(static_cast<std::uint16_t>(loc2),
                                         static_cast<std::uint16_t>(loc1)));
      }
      penalty += cost;
    }
  }

  return penalty;
}
