//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Definitions.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/OpType.hpp"
#include "sc/Architecture.hpp"
#include "sc/DataLogger.hpp"
#include "sc/configuration/AvailableArchitecture.hpp"
#include "sc/configuration/EarlyTermination.hpp"
#include "sc/configuration/Heuristic.hpp"
#include "sc/configuration/InitialLayout.hpp"
#include "sc/configuration/Layering.hpp"
#include "sc/configuration/LookaheadHeuristic.hpp"
#include "sc/configuration/Method.hpp"
#include "sc/heuristic/HeuristicMapper.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <iterator>
#include <memory>
#include <nlohmann/json.hpp>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

constexpr qc::OpType SWAP = qc::OpType::SWAP;
constexpr double FLOAT_TOLERANCE = 1e-6;

/**
 * @brief Get id of the final node in a given layer from a data log.
 */
std::size_t getFinalNodeFromDatalog(std::string dataLoggingPath,
                                    std::size_t layer) {
  if (dataLoggingPath.back() != '/') {
    dataLoggingPath += '/';
  }
  auto layerFile = std::ifstream(dataLoggingPath + "layer_" +
                                 std::to_string(layer) + ".json");
  if (!layerFile.is_open()) {
    throw std::runtime_error("Could not open file " + dataLoggingPath +
                             "layer_" + std::to_string(layer) + ".json");
  }
  const auto layerJson = nlohmann::basic_json<>::parse(layerFile);
  if (layerJson.find("final_node_id") == layerJson.end()) {
    throw std::runtime_error("Missing key \"final_node_id\" in " +
                             dataLoggingPath + "layer_" +
                             std::to_string(layer) + ".json");
  }
  return layerJson["final_node_id"];
}

/**
 * @brief parses all nodes in a given layer from a data log and enter them
 * into `nodes` with each node at the position corresponding to its id.
 *
 * Only logged values are entered into the nodes, all other values are left at
 * default (e.g. `validMappedTwoQubitGates` and `sharedSwaps`)
 */
void parseNodesFromDatalog(std::string dataLoggingPath, std::size_t layer,
                           std::vector<HeuristicMapper::Node>& nodes) {
  if (dataLoggingPath.back() != '/') {
    dataLoggingPath += '/';
  }
  const std::string layerNodeFilePath =
      dataLoggingPath + "nodes_layer_" + std::to_string(layer) + ".csv";
  auto layerNodeFile = std::ifstream(layerNodeFilePath);
  if (!layerNodeFile.is_open()) {
    throw std::runtime_error("Could not open file " + layerNodeFilePath);
  }
  // iterating over the lines in the csv file and then over the entries
  // separated by ';'
  std::string line;
  while (std::getline(layerNodeFile, line)) {
    if (line.empty()) {
      continue;
    }
    std::string col;
    std::size_t nodeId = 0;

    std::stringstream lineStream(line);
    if (std::getline(lineStream, col, ';')) {
      nodeId = std::stoull(col);
      if (nodeId >= nodes.size()) {
        throw std::runtime_error("Node id " + std::to_string(nodeId) +
                                 " out of range in " + layerNodeFilePath);
      }
      nodes[nodeId].id = nodeId;
    } else {
      throw std::runtime_error("Missing value for node id in " +
                               layerNodeFilePath);
    }
    auto& node = nodes[nodeId];
    if (std::getline(lineStream, col, ';')) {
      node.parent = std::stoull(col);
    } else {
      throw std::runtime_error("Missing value for parent node id in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      node.costFixed = std::stod(col);
    } else {
      throw std::runtime_error("Missing value for fixed cost in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      node.costHeur = std::stod(col);
    } else {
      throw std::runtime_error("Missing value for heuristic cost in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      node.lookaheadPenalty = std::stod(col);
    } else {
      throw std::runtime_error("Missing value for lookahead penalty in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      const std::size_t validMapping = std::stoull(col);
      if (validMapping > 1) {
        throw std::runtime_error("Non-boolean value " +
                                 std::to_string(validMapping) +
                                 " for validMapping in " + layerNodeFilePath);
      }
      node.validMapping = static_cast<bool>(validMapping);
    } else {
      throw std::runtime_error("Missing value for validMapping in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      node.depth = std::stoull(col);
    } else {
      throw std::runtime_error("Missing value for depth in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      std::stringstream qubitMapBuffer(col);
      std::string entry;
      for (std::size_t i = 0; std::getline(qubitMapBuffer, entry, ','); ++i) {
        auto qubit = static_cast<std::int16_t>(std::stoi(entry));
        node.qubits.at(i) = qubit;
        if (qubit >= 0) {
          node.locations.at(static_cast<std::size_t>(qubit)) =
              static_cast<std::int16_t>(i);
        }
      }
    } else {
      throw std::runtime_error("Missing value for qubit layout in " +
                               layerNodeFilePath);
    }
    if (std::getline(lineStream, col, ';')) {
      std::stringstream swapBuffer(col);
      std::string entry;
      while (std::getline(swapBuffer, entry, ',')) {
        std::uint16_t q1 = 0;
        std::uint16_t q2 = 0;
        std::string opTypeStr;
        std::stringstream(entry) >> q1 >> q2 >> opTypeStr;
        qc::OpType opType = SWAP;
        if (!opTypeStr.empty()) {
          // if no opType is given, the default value is SWAP
          opType = qc::opTypeFromString(opTypeStr);
        }
        node.swaps.emplace_back(q1, q2, opType);
      }
    }
  }
}

/**
 * @brief Get the path from a node to the root node (id of the given node is the
 * first element, id of the root is last)
 *
 * @param nodes vector of all nodes (each at the position corresponding to its
 * id)
 * @param nodeId id of the node from which to find the path to the root
 */
std::vector<std::size_t>
getPathToRoot(std::vector<HeuristicMapper::Node>& nodes, std::size_t nodeId) {
  std::vector<std::size_t> path{};
  if (nodeId >= nodes.size() || nodes[nodeId].id != nodeId) {
    throw std::runtime_error("Invalid node id " + std::to_string(nodeId));
  }
  auto* node = &nodes[nodeId];
  while (node->parent != node->id) {
    path.push_back(node->id);
    if (node->parent >= nodes.size() ||
        nodes[node->parent].id != node->parent) {
      throw std::runtime_error("Invalid parent id " +
                               std::to_string(node->parent) + " for node " +
                               std::to_string(node->id));
    }
    node = &nodes[node->parent];
  }
  path.push_back(node->id);
  return path;
}

class InternalsTest : public HeuristicMapper, public testing::Test {
protected:
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static Architecture defaultArch;

  InternalsTest() : HeuristicMapper(qc::QuantumComputation{1}, defaultArch) {}
  void SetUp() override { results = MappingResults{}; }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
Architecture InternalsTest::defaultArch{1, {}};

TEST_F(InternalsTest, NodeCostCalculation) {
  results.config.heuristic = Heuristic::GateCountMaxDistance;
  results.config.lookaheadHeuristic = LookaheadHeuristic::None;
  results.config.layering = Layering::Disjoint2qBlocks;

  architecture->loadCouplingMap(5, {{0, 1}, {1, 2}, {3, 1}, {4, 3}});
  qc = qc::QuantumComputation{5};
  // layer 0 distances:
  // 0-1: 2 swaps & 2 reversals
  // 2-3: 0 swaps & 1 reversals
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(1, 0);
  qc.cx(1, 0);
  qc.cx(3, 2);
  // Architecture::distance currently counts at most 1 reversal per qubit pair
  createLayers();

  EXPECT_EQ(layers.size(), 1)
      << "layering failed, not able to test node cost calculation";

  const std::vector<Exchange> swaps{Exchange(0, 1, qc::OpType::Teleportation),
                                    Exchange(1, 2, SWAP)};

  HeuristicMapper::Node node(0, 0, {4, 3, 1, 2, 0}, {4, 2, 3, 1, 0}, swaps,
                             {{2, 3}}, 5., 0);
  EXPECT_NEAR(node.costFixed, 5., FLOAT_TOLERANCE);
  EXPECT_NEAR(node.lookaheadPenalty, 0., FLOAT_TOLERANCE);
  EXPECT_EQ(node.validMappedTwoQubitGates.size(), 1);

  results.config.heuristic = Heuristic::GateCountSumDistance;
  updateHeuristicCost(0, node);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE * 2,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.costFixed, 5., FLOAT_TOLERANCE)
      << "updateHeuristicCost should not change costFixed";

  results.config.heuristic = Heuristic::GateCountMaxDistance;
  updateHeuristicCost(0, node);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.costFixed, 5., FLOAT_TOLERANCE)
      << "updateHeuristicCost should not change costFixed";

  applySWAP({3, 4}, 0, node);
  EXPECT_NEAR(node.costFixed, 5. + COST_UNIDIRECTIONAL_SWAP, FLOAT_TOLERANCE);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_EQ(node.validMappedTwoQubitGates.size(), 0);

  node.lookaheadPenalty = 0.;
  EXPECT_NEAR(node.getTotalCost(),
              5. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.getTotalFixedCost(), 5. + COST_UNIDIRECTIONAL_SWAP,
              FLOAT_TOLERANCE);

  node.lookaheadPenalty = 2.;
  EXPECT_NEAR(node.getTotalCost(),
              7. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.getTotalFixedCost(), 7. + COST_UNIDIRECTIONAL_SWAP,
              FLOAT_TOLERANCE);

  recalculateFixedCost(0, node);
  EXPECT_EQ(node.validMappedTwoQubitGates.size(), 0);
  EXPECT_NEAR(node.costFixed, COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.getTotalCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 3 +
                  COST_DIRECTION_REVERSE,
              FLOAT_TOLERANCE);
  EXPECT_NEAR(node.getTotalFixedCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              FLOAT_TOLERANCE);
}

TEST_F(InternalsTest, NodeLookaheadCalculation) {
  results.config.heuristic = Heuristic::GateCountMaxDistance;
  results.config.lookaheadHeuristic = LookaheadHeuristic::None;
  results.config.layering = Layering::Disjoint2qBlocks;

  architecture->loadCouplingMap(5, {{0, 1}, {1, 2}, {3, 1}, {4, 3}});
  qc = qc::QuantumComputation{5};
  // layer 0 distances:
  // 0-1: 2 swaps & 2 reversals
  // 2-3: 0 swaps & 1 reversals
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(0, 1);
  qc.cx(1, 0);
  qc.cx(1, 0);
  qc.cx(3, 2);

  // layer 1 distances:
  // 0-4: 2 swaps & 0 reversals
  // 1-2: 1 swaps & 1 reversals
  qc.cx(0, 4);
  qc.cx(0, 4);
  qc.cx(1, 2);

  // layer 2 distances:
  // 0-1: 2 swaps & 0 reversals
  qc.cx(0, 1);

  // layer 3 distances:
  // 0-4: 2 swaps & 0 reversals
  qc.cx(4, 0);

  // Architecture::distance currently counts at most 1 reversal per qubit pair

  createLayers();

  EXPECT_EQ(layers.size(), 4)
      << "layering failed, not able to test node cost calculation";

  const std::vector<Exchange> swaps{Exchange(0, 1, qc::OpType::Teleportation),
                                    Exchange(1, 2, SWAP)};

  HeuristicMapper::Node node(0, 0, {4, 3, 1, 2, 0}, {4, 2, 3, 1, 0}, swaps,
                             {{2, 3}}, 5., 0);
  EXPECT_NEAR(node.lookaheadPenalty, 0., FLOAT_TOLERANCE);

  results.config.firstLookaheadFactor = 0.75;
  results.config.lookaheadFactor = 0.5;

  results.config.lookaheadHeuristic = LookaheadHeuristic::GateCountMaxDistance;
  results.config.nrLookaheads = 1;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty, 0.75 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 2;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 3;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 4;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);

  results.config.lookaheadHeuristic = LookaheadHeuristic::GateCountSumDistance;
  results.config.nrLookaheads = 1;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (3 * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 2;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (3 * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 3;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (3 * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);
  results.config.nrLookaheads = 4;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (3 * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE) +
                  0.75 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP) +
                  0.75 * 0.5 * 0.5 * (2 * COST_UNIDIRECTIONAL_SWAP),
              FLOAT_TOLERANCE);

  node.qubits = {4, 3, 1, -1, -1};
  node.locations = {-1, 2, -1, 1, 0};

  results.config.lookaheadHeuristic = LookaheadHeuristic::GateCountMaxDistance;
  results.config.nrLookaheads = 1;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE),
              FLOAT_TOLERANCE);

  results.config.lookaheadHeuristic = LookaheadHeuristic::GateCountSumDistance;
  results.config.nrLookaheads = 1;
  updateLookaheadPenalty(0, node);
  EXPECT_NEAR(node.lookaheadPenalty,
              0.75 * (2 * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE),
              FLOAT_TOLERANCE);
}

class TestHeuristics
    : public testing::TestWithParam<std::tuple<Heuristic, std::string>> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";
  std::string testCalibrationDir = "../../../extern/calibration/";

  qc::QuantumComputation qc;
  std::string circuitName;
  Architecture ibmqYorktown; // 5 qubits
  Architecture ibmqLondon;   // 5 qubits (with calibration)
  std::unique_ptr<HeuristicMapper> ibmqYorktownMapper;
  std::unique_ptr<HeuristicMapper> ibmqLondonMapper;
  Architecture ibmQX5; // 16 qubits
  std::unique_ptr<HeuristicMapper> ibmQX5Mapper;
  Configuration settings{};

  static const std::unordered_map<std::string,
                                  std::vector<std::vector<std::int16_t>>>
      OPTIMAL_SOLUTIONS;

  void SetUp() override {
    std::string cn = std::get<1>(GetParam());
    std::replace(cn.begin(), cn.end(), '-', '_');
    std::stringstream ss{};
    ss << cn << "_" << toString(std::get<0>(GetParam()));
    const std::string testName = ss.str();

    circuitName = std::get<1>(GetParam());
    qc.import(testExampleDir + circuitName + ".qasm");
    ibmqYorktown.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
    ibmqYorktownMapper = std::make_unique<HeuristicMapper>(qc, ibmqYorktown);
    ibmqLondonMapper = std::make_unique<HeuristicMapper>(qc, ibmqLondon);
    ibmQX5.loadCouplingMap(AvailableArchitecture::IbmQx5);
    ibmQX5Mapper = std::make_unique<HeuristicMapper>(qc, ibmQX5);
    settings.debug = true;
    settings.automaticLayerSplits = false;
    settings.initialLayout = InitialLayout::Identity;
    settings.layering = Layering::Disjoint2qBlocks;
    settings.lookaheadHeuristic = LookaheadHeuristic::None;
    settings.heuristic = std::get<0>(GetParam());
    settings.dataLoggingPath = "test_log/heur_properties_" + testName + "_" +
                               std::to_string(std::random_device{}());
  }

  void TearDown() override {
    std::filesystem::remove_all(settings.dataLoggingPath);
  }
};

const std::unordered_map<std::string, std::vector<std::vector<std::int16_t>>>
    TestHeuristics::OPTIMAL_SOLUTIONS{
        {"3_17_13",
         {{0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {1, 0, 2},
          {1, 2, 0},
          {1, 2, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 0, 1},
          {2, 1, 0},
          {2, 0, 1},
          {2, 0, 1},
          {2, 1, 0},
          {1, 0, 2},
          {1, 2, 0},
          {1, 2, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {2, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {0, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2},
          {0, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2},
          {0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2},
          {2, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0}}},
        {"ex-1_166",
         {{0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {1, 0, 2},
          {1, 0, 2},
          {1, 2, 0},
          {1, 2, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {0, 1, 2},
          {0, 1, 2},
          {1, 0, 2},
          {1, 0, 2},
          {0, 1, 2},
          {1, 0, 2},
          {1, 0, 2},
          {1, 2, 0},
          {1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0}}},
        {"ham3_102",
         {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2},
          {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {1, 0, 2}, {1, 0, 2},
          {1, 2, 0}, {1, 2, 0}, {2, 1, 0}, {2, 1, 0}, {2, 0, 1}, {2, 1, 0},
          {0, 1, 2}, {1, 0, 2}, {1, 0, 2}, {0, 1, 2}, {1, 0, 2}, {1, 0, 2},
          {1, 2, 0}, {1, 2, 0}, {1, 2, 0}}},
        {"miller_11",
         {{0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {0, 1, 2},
          {1, 0, 2},
          {1, 2, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {2, 1, 0},
          {0, 1, 2},
          {0, 1, 2},
          {0, 2, 1},
          {0, 2, 1},
          {0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1},
          {0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1},
          {2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1},
          {2, 0, 1},
          {2, 1, 0},
          {2, 1, 0},
          {1, 2, 0},
          {1, 2, 0},
          {1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {1, 2, 0},
          {1, 2, 0},
          {1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {1, 2, 0},
          {1, 2, 0},
          {2, 1, 0},
          {1, 2, 0},
          {1, 2, 0},
          {1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0},
          {1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0}}},
        {"4gt11_84", {{0, 1, 2, -1, 4}, {0, 1, 2, -1, 4}, {0, 1, 2, -1, 4},
                      {0, 1, 2, -1, 4}, {0, 1, 2, -1, 4}, {0, 1, 2, -1, 4},
                      {2, 1, 0, 3, 4},  {0, 1, 2, -1, 4}, {0, 1, 2, -1, 4},
                      {1, 0, 2, -1, 4}, {1, 2, 0, -1, 4}, {2, 1, 0, -1, 4},
                      {2, 0, 1, -1, 4}, {4, 2, 1, 0, 3},  {0, 1, 2, -1, 4},
                      {0, 1, 2, -1, 4}, {0, 2, 1, -1, 4}, {0, 2, 1, -1, 4},
                      {2, 0, 1, -1, 4}, {2, 0, 1, -1, 4}, {2, 1, 0, 3, 4}}},
        {"4mod5-v0_20",
         {{0, 2, 1, 3, 4}, {0, 2, 1, 3, 4}, {0, 2, 1, 3, 4}, {0, 2, 4, 3, 1},
          {0, 4, 2, 3, 1}, {0, 4, 1, 3, 2}, {0, 4, 2, 3, 1}, {0, 4, 2, 3, 1},
          {4, 0, 2, 1, 3}, {4, 2, 0, 1, 3}, {4, 1, 0, 2, 3}, {4, 2, 0, 1, 3},
          {4, 2, 0, 1, 3}, {4, 1, 0, 2, 3}, {4, 2, 0, 1, 3}, {4, 2, 0, 1, 3},
          {0, 2, 1, 3, 4}, {0, 2, 1, 3, 4}, {0, 2, 3, 1, 4}, {0, 3, 2, 4, 1},
          {0, 3, 4, 2, 1}, {0, 3, 4, 1, 2}, {0, 3, 4, 2, 1}, {0, 3, 4, 2, 1}}},
        {"mod5d1_63",
         {{0, 2, 1, 3, 4},
          {0, 2, 1, 3, 4},
          {1, 2, 0, 3, 4},
          {1, 2, 0, 3, 4},
          {1, 2, 0, 3, 4},
          {1, 2, 4, 3, 0},
          {4, 2, 1, 3, 0},
          {4, 2, 0, 3, 1},
          {4, 2, 1, 3, 0},
          {4, 2, 1, 3, 0},
          {4, 2, 0, 3, 1},
          {4, 0, 2, 1, 3},
          {4, 1, 2, 0, 3},
          {4, 0, 2, 1, 3},
          {4, 0, 2, 1, 3},
          {4, 0, 2, 1, 3},
          {4, 1, 2, 0, 3},
          {4, 1, 2, 0, 3},
          {4, 0, 2, 1, 3},
          {4, 1, 2, 0, 3},
          {4, 1, 2, 0, 3},
          {4, 0, 2, 1, 3},
          {0, 2, 1, 3, 4},
          {0, 2, 3, 1, 4},
          {2, 3, 1, 0, 4},
          {2, 3, 1, 0, 4},
          {2, 3, 1, 0, 4},
          {2, 3, 1, 4, 0},
          {2, 3, 4, 1, 0},
          {2, 3, 4, 0, 1},
          {2, 3, 4, 1, 0},
          {2, 3, 4, 1, 0},
          {-1, 3, 1, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4}}},
        {"ising_model_10",
         {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}}},
        {"rd73_140",
         {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
          {0, 2, 3, 4, 5, 1, 7, 6, 8, 9},
          {2, 3, 4, 0, 7, 5, 1, 6, 8, 9},
          {2, 3, 4, 0, 1, 7, 5, 6, 8, 9},
          {2, 3, 4, 0, 1, 7, 5, 6, 8, 9},
          {2, 3, 4, 0, 7, 1, 5, 6, 8, 9},
          {2, 7, 3, 4, 0, 1, 5, 6, 8, 9},
          {2, 3, 4, 0, 7, 8, 1, 5, 6, 9},
          {3, 4, 2, 8, 0, 7, 1, 5, 6, 9},
          {3, 4, 8, 2, 7, 0, 1, 5, 6, 9},
          {3, 4, 8, 7, 2, 0, 1, 5, 6, 9},
          {3, 4, 8, 2, 7, 0, 1, 5, 6, 9},
          {3, 8, 4, 2, 7, 0, 1, 5, 6, 9},
          {3, 8, 9, 4, 7, 2, 1, 5, 6, -1, -1, -1, 0},
          {8, 3, 9, 4, 7, 2, 1, 5, 6, -1, -1, -1, 0},
          {8, 3, 9, 4, 7, 1, 2, 5, 6, -1, -1, -1, 0},
          {8, 9, 3, 4, 7, 1, 2, 5, 6, -1, -1, -1, 0},
          {8, 9, 3, 4, 7, 2, 1, 5, 6, -1, -1, -1, 0},
          {9, 8, 3, 4, 7, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 4, 3, 7, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 7, 4, 3, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 3, 7, 4, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 3, 7, 4, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 7, 3, 4, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 3, 7, 4, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 4, 3, 7, 1, 2, 5, 6, -1, -1, -1, 0},
          {9, 8, 4, 3, 2, 7, 1, 5, 6, -1, -1, -1, 0},
          {9, 8, -1, 2, 3, 7, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {9, -1, 8, 2, 7, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {4, -1, 8, 2, 3, 7, 1, 5, 6, -1, -1, -1, 0, -1, -1, 9},
          {4, -1, 8, 2, 3, 7, 1, 5, 6, -1, -1, -1, 0, -1, -1, 9},
          {4, 8, -1, 2, 7, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 9},
          {9, 8, -1, 7, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, 4},
          {9, 8, 7, -1, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, 4},
          {9, 7, 8, -1, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {9, 8, 7, -1, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {9, 8, 7, -1, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {9, 7, 8, -1, 2, 3, 1, 5, 6, -1, -1, -1, 0, -1, -1, 4},
          {4, 7, 8, 5, -1, 2, 3, 1, 6, -1, -1, -1, 0, -1, -1, 9},
          {9, 7, 5, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 8, 7, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 4},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 4, 5, 7, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 8, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 7},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 7, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 8},
          {9, 5, 8, 4, 3, -1, 2, 1, 6, -1, -1, -1, 0, -1, -1, 7},
          {5, 9, 8, 6, 4, 3, -1, 2, 1, -1, -1, -1, 0, -1, -1, 7},
          {7, 9, 8, 6, 4, 3, -1, 2, 1, -1, -1, -1, 0, 5},
          {-1, 9, 6, 8, 4, 3, -1, 2, 1, -1, -1, -1, 0, -1, 5, 7},
          {-1, 9, 6, 8, -1, 3, -1, 2, 1, -1, -1, -1, 0, 4, 7, 5},
          {-1, 9, 8, 6, -1, 3, -1, 2, 1, -1, -1, -1, 0, 4, 5, 7},
          {-1, 9, 6, 8, -1, 3, -1, 2, 1, -1, -1, -1, 0, 4, 5, 7},
          {-1, 9, 6, 8, -1, 3, -1, 2, 1, -1, -1, -1, 0, 4, 7, 5},
          {-1, 9, 6, 8, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 4, 5},
          {-1, 9, 6, 8, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 4, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 6, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 6, 7, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 6, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 6, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 6, 7, 5},
          {-1, 9, 8, 4, -1, 3, -1, 2, 1, -1, -1, -1, 0, 7, 6, 5}}}};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, TestHeuristics,
    testing::Combine(
        testing::Values(
            Heuristic::GateCountMaxDistance, Heuristic::GateCountSumDistance,
            Heuristic::GateCountSumDistanceMinusSharedSwaps,
            Heuristic::GateCountMaxDistanceOrSumDistanceMinusSharedSwaps,
            Heuristic::FidelityBestLocation),
        testing::Values("3_17_13",        // 5q
                        "ex-1_166",       // 5q
                        "ham3_102",       // 5q
                        "miller_11",      // 5q
                        "4gt11_84",       // 5q
                        "4mod5-v0_20",    // 5q
                        "mod5d1_63",      // 5q
                        "ising_model_10", // 16q
                        "rd73_140"        // 16q
                        )),
    [](const testing::TestParamInfo<TestHeuristics::ParamType>& inf) {
      std::string name = std::get<1>(inf.param);
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name << "_" << toString(std::get<0>(inf.param));
      return ss.str();
    });

TEST_P(TestHeuristics, HeuristicProperties) {
  EXPECT_TRUE(!isAdmissible(settings.heuristic) ||
              isPrincipallyAdmissible(settings.heuristic))
      << "Admissible heuristics are by definition also principally admissible: "
      << toString(settings.heuristic);

  EXPECT_TRUE(
      !(isTight(settings.heuristic) && isFidelityAware(settings.heuristic)))
      << "Fidelity-aware heuristics cannot be tight because of the "
         "non-convexity of the fidelity-aware cost function: "
      << toString(settings.heuristic);

  // list of nodes for each search process (i.e. each layer) in all mappings
  // each node is at the position corresponding to its id; positions of unused
  // ids are filled with default values (i.e. node.id = 0)
  std::vector<std::vector<HeuristicMapper::Node>> allNodes{};
  std::vector<std::string> layerNames{};
  std::vector<std::size_t> finalSolutionIds{};

  // map to IBM Yorktown if possible
  if (qc.getNqubits() <= ibmqYorktown.getNqubits()) {
    if (isFidelityAware(settings.heuristic)) {
      EXPECT_THROW(ibmqYorktownMapper->map(settings), QMAPException);
    } else {
      ibmqYorktownMapper->map(settings);
      auto results = ibmqYorktownMapper->getResults();
      for (std::size_t i = 0; i < results.layerHeuristicBenchmark.size(); ++i) {
        auto& nodes = allNodes.emplace_back(
            results.layerHeuristicBenchmark.at(i).generatedNodes,
            HeuristicMapper::Node{ibmqYorktown.getNqubits(), 0});
        layerNames.emplace_back("on ibmq_yorktown in layer " +
                                std::to_string(i));
        parseNodesFromDatalog(settings.dataLoggingPath, i, nodes);
        finalSolutionIds.emplace_back(
            getFinalNodeFromDatalog(settings.dataLoggingPath, i));
      }
    }
  }

  // map to IBM London if possible
  if (qc.getNqubits() <= ibmqLondon.getNqubits()) {
    ibmqLondonMapper->map(settings);
    auto results = ibmqLondonMapper->getResults();
    for (std::size_t i = 0; i < results.layerHeuristicBenchmark.size(); ++i) {
      auto& nodes = allNodes.emplace_back(
          results.layerHeuristicBenchmark.at(i).generatedNodes,
          HeuristicMapper::Node{ibmqLondon.getNqubits(), 0});
      layerNames.emplace_back("on ibmq_london in layer " + std::to_string(i));
      parseNodesFromDatalog(settings.dataLoggingPath, i, nodes);
      finalSolutionIds.emplace_back(
          getFinalNodeFromDatalog(settings.dataLoggingPath, i));
    }
  }

  // map to IBM QX5 if possible
  if (qc.getNqubits() <= ibmQX5.getNqubits()) {
    if (isFidelityAware(settings.heuristic)) {
      EXPECT_THROW(ibmQX5Mapper->map(settings), QMAPException);
    } else {
      ibmQX5Mapper->map(settings);
      auto results = ibmQX5Mapper->getResults();
      for (std::size_t i = 0; i < results.layerHeuristicBenchmark.size(); ++i) {
        auto& nodes = allNodes.emplace_back(
            results.layerHeuristicBenchmark.at(i).generatedNodes,
            HeuristicMapper::Node{ibmQX5.getNqubits(), 0});
        layerNames.emplace_back("on ibmQX5 in layer " + std::to_string(i));
        parseNodesFromDatalog(settings.dataLoggingPath, i, nodes);
        finalSolutionIds.emplace_back(
            getFinalNodeFromDatalog(settings.dataLoggingPath, i));
      }
    }
  }

  for (std::size_t i = 0; i < allNodes.size(); ++i) {
    auto& nodes = allNodes.at(i);
    auto& finalSolutionId = finalSolutionIds.at(i);

    if (finalSolutionId >= nodes.size()) {
      FAIL() << "Final solution id " << finalSolutionId << " out of bounds "
             << layerNames.at(i) << " (number of nodes: " << nodes.size()
             << ")";
    }
    if (nodes.at(finalSolutionId).id != finalSolutionId) {
      FAIL() << "Node id " << nodes.at(finalSolutionId).id
             << " does not match final solution id " << finalSolutionId << " "
             << layerNames.at(i);
    }
    auto& finalSolutionNode = nodes.at(finalSolutionId);
    EXPECT_TRUE(finalSolutionNode.validMapping);

    if (isPrincipallyAdmissible(settings.heuristic)) {
      // for principally admissible heuristics all nodes on the optimal
      // solution path should have
      // node.costFixed+node.costHeur <= finalSolutionNode.costFixed
      auto solutionPath = getPathToRoot(nodes, finalSolutionId);
      for (auto nodeId : solutionPath) {
        if (nodes.at(nodeId).id != nodeId) {
          throw std::runtime_error("Invalid node id " + std::to_string(nodeId) +
                                   " " + layerNames.at(i));
        }
        EXPECT_LE(nodes.at(nodeId).getTotalCost(), finalSolutionNode.costFixed)
            << "Heuristic " << toString(settings.heuristic)
            << " is not principally admissible " << layerNames.at(i)
            << " in node " << nodeId;
      }
      if (isTight(settings.heuristic)) {
        // Principally admissible heuristics are guaranteed to find a solution
        // with optimal cost (given lookahead is disabled). If there are
        // multiple optimal solutions, though, and the heuristic is not tight,
        // the differences in heuristics can lead to a different node ordering
        // and thereby different solutions (potentially even resulting in
        // different costs in later layers)
        // However, if a heuristic is both principally admissible and tight,
        // it is guaranteed to always find the same solution as any other such
        // heuristic.
        if (OPTIMAL_SOLUTIONS.find(circuitName) == OPTIMAL_SOLUTIONS.end() ||
            OPTIMAL_SOLUTIONS.at(circuitName).size() <= i) {
          throw std::runtime_error(
              "Missing precalculated optimal solution for circuit " +
              circuitName);
        }
        std::size_t finalLayoutLastIndex = 0;
        for (std::size_t j = finalSolutionNode.qubits.size() - 1; j > 0; --j) {
          if (finalSolutionNode.qubits.at(j) != -1) {
            finalLayoutLastIndex = j;
            break;
          }
        }
        std::vector<std::int16_t> finalLayout{};
        std::copy(finalSolutionNode.qubits.begin(),
                  finalSolutionNode.qubits.begin() +
                      static_cast<std::vector<int16_t>::difference_type>(
                          finalLayoutLastIndex + 1),
                  std::back_inserter(finalLayout));
        EXPECT_EQ(finalLayout, OPTIMAL_SOLUTIONS.at(circuitName).at(i))
            << "Heuristic " << toString(settings.heuristic)
            << " did not find the optimal solution " << layerNames.at(i);
      }
    }

    for (std::size_t j = 0; j < nodes.size(); ++j) {
      const auto& node = nodes.at(j);
      if (j != node.id) {
        continue;
      }

      if (isNonDecreasing(settings.heuristic)) {
        if (node.parent != node.id) {
          if (node.parent >= nodes.size() ||
              nodes.at(node.parent).id != node.parent) {
            FAIL() << "Invalid parent id " << node.parent << " for node "
                   << node.id << " " << layerNames.at(i);
          }
          EXPECT_GE(node.getTotalCost(), nodes.at(node.parent).getTotalCost())
              << "Heuristic " << toString(settings.heuristic)
              << " does not result in non-decreasing cost estimation "
              << layerNames.at(i) << " in node " << node.id;
        }
      }

      EXPECT_NEAR(node.lookaheadPenalty, 0., FLOAT_TOLERANCE)
          << "Lookahead penalty not 0 " << layerNames.at(i)
          << " even though lookahead has been deactivated";

      if (node.validMapping) {
        if (isTight(settings.heuristic)) {
          // tight heuristics are 0 in any goal node
          EXPECT_NEAR(node.costHeur, 0., FLOAT_TOLERANCE)
              << "Heuristic " << toString(settings.heuristic)
              << " is not tight " << layerNames.at(i) << " in node " << node.id;
        }

        if (isAdmissible(settings.heuristic)) {
          // for admissible heuristics all nodes should have
          // node.costFixed+node.costHeur <= solutionNode.costFixed,
          // where solutionNode is the best goal node reachable from the node;
          // since reachability in a directed tree is equivalent to the being
          // on the same path to the root, one can also check that all nodes on
          // the path to the root from any goal node fulfill this condition
          auto path = getPathToRoot(nodes, node.id);
          for (auto nodeId : path) {
            auto& n = nodes.at(nodeId);
            EXPECT_LE(n.getTotalCost(), node.costFixed)
                << "Heuristic " << toString(settings.heuristic)
                << " is not admissible " << layerNames.at(i) << " in node "
                << nodeId;
          }
        }
      }
    }
  }
}

TEST(Functionality, HeuristicBenchmark) {
  /*
      3
     / \
    4   2
    |   |
    0---1
  */
  Architecture architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {3, 4}, {4, 3}, {4, 0}, {0, 4}};
  architecture.loadCouplingMap(5, cm);

  qc::QuantumComputation qc{5, 5};
  qc.cx(4, 2);
  qc.cx(3, 1);
  qc.cx(4, 1);

  qc.barrier({0, 1, 2, 3, 4});
  for (size_t i = 0; i < 5; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  const auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);
  Configuration settings{};
  settings.heuristic = Heuristic::GateCountMaxDistance;
  settings.layering = Layering::DisjointQubits;
  settings.initialLayout = InitialLayout::Identity;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.debug = true;
  mapper->map(settings);
  auto& result = mapper->getResults();

  /*
  generated nodes (unit of costs: COST_BIDIRECTIONAL_SWAP):
  layer 1:
    0: {swaps: {}, cost: 0, heur: 1, total: 1}
  --- priority queue: [0] -> expand: 0
    1: {swaps: {{0, 1}}, cost: 1, heur: 1, total: 2}
    2: {swaps: {{1, 2}}, cost: 1, heur: 1, total: 2}
    3: {swaps: {{2, 3}}, cost: 1, heur: 0, total: 1}
    4: {swaps: {{3, 4}}, cost: 1, heur: 1, total: 2}
    5: {swaps: {{4, 0}}, cost: 1, heur: 1, total: 2}
  --- priority queue: [3,1,2,4,5] -> done: 3
  layer 2:
    0: {swaps: {}, cost: 0, heur: 1, total: 1}
  --- priority queue: [0] -> expand: 0
    1: {swaps: {{0, 1}}, cost: 1, heur: 0, total: 1}
    2: {swaps: {{1, 2}}, cost: 1, heur: 1, total: 2}
    3: {swaps: {{3, 4}}, cost: 1, heur: 1, total: 2}
    4: {swaps: {{4, 0}}, cost: 1, heur: 0, total: 1}
  --- priority queue: [1,4,2,3] -> done: 1
  */

  EXPECT_EQ(result.layerHeuristicBenchmark.size(), 2);
  const auto& layerResults0 = result.layerHeuristicBenchmark[0];
  EXPECT_EQ(layerResults0.solutionDepth, 1);
  EXPECT_EQ(layerResults0.generatedNodes, 6);
  EXPECT_EQ(layerResults0.expandedNodes, 1);
  EXPECT_NEAR(layerResults0.averageBranchingFactor, 5.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(layerResults0.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  const auto& layerResults1 = result.layerHeuristicBenchmark[1];
  EXPECT_EQ(layerResults1.solutionDepth, 1);
  EXPECT_EQ(layerResults1.generatedNodes, 5);
  EXPECT_EQ(layerResults1.expandedNodes, 1);
  EXPECT_NEAR(layerResults1.averageBranchingFactor, 4.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(layerResults1.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_EQ(result.heuristicBenchmark.generatedNodes, 11);
  EXPECT_EQ(result.heuristicBenchmark.expandedNodes, 2);
  EXPECT_NEAR(result.heuristicBenchmark.averageBranchingFactor, 4.5,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(result.heuristicBenchmark.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
}

TEST(Functionality, EmptyDump) {
  qc::QuantumComputation qc{1};
  qc.x(0);
  Architecture arch{1, {}};
  HeuristicMapper mapper(qc, arch);
  mapper.dumpResult("test.qasm");
  mapper.map({});
  EXPECT_NO_THROW(mapper.dumpResult("test.qasm"););
  EXPECT_THROW(mapper.dumpResult("test.dummy"), QMAPException);
  std::filesystem::remove("test.qasm");
}

TEST(Functionality, BenchmarkGeneratedNodes) {
  qc::QuantumComputation qc{16, 16};
  qc.cx(0, 6);
  for (std::size_t i = 0; i < 16; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }
  Architecture ibmQX5{};
  ibmQX5.loadCouplingMap(AvailableArchitecture::IbmQx5);
  auto ibmQX5Mapper = std::make_unique<HeuristicMapper>(qc, ibmQX5);

  Configuration settings{};
  settings.heuristic = Heuristic::GateCountMaxDistance;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.layering = Layering::IndividualGates;
  settings.automaticLayerSplits = false;
  settings.initialLayout = InitialLayout::Identity;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.useTeleportation = false;
  settings.debug = true;
  ibmQX5Mapper->map(settings);
  auto results = ibmQX5Mapper->getResults();

  EXPECT_EQ(results.heuristicBenchmark.generatedNodes, 30);
  EXPECT_EQ(results.layerHeuristicBenchmark.at(0).generatedNodes, 30);
}

TEST(Functionality, InvalidSettings) {
  qc::QuantumComputation qc{1};
  qc.x(0);
  Architecture arch{1, {}};
  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.1);
  arch.loadProperties(props);
  HeuristicMapper mapper(qc, arch);
  auto config = Configuration{};
  config.method = Method::Heuristic;
  config.heuristic = Heuristic::GateCountMaxDistance;
  config.lookaheadHeuristic = LookaheadHeuristic::GateCountMaxDistance;
  // invalid layering
  config.layering = Layering::OddGates;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.layering = Layering::QubitTriangle;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.layering = Layering::IndividualGates;
  // fidelity-aware heuristic with non-fidelity-aware lookahead heuristic
  config.heuristic = Heuristic::FidelityBestLocation;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.lookaheadHeuristic = LookaheadHeuristic::None;
  // fidelity-aware heuristic with teleportation
  config.teleportationQubits = 2;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.teleportationQubits = 0;
  // valid settings
  EXPECT_NO_THROW(mapper.map(config));
}

TEST(Functionality, NoMeasurementsAdded) {
  using namespace qc::literals;
  // construct circuit
  qc::QuantumComputation qc{4U};
  qc.cx(0_pc, 1);
  qc.cx(2_pc, 1);
  qc.cx(3_pc, 1);

  // load architecture
  Architecture arch{};
  arch.loadCouplingMap(AvailableArchitecture::IbmqLondon);

  // create heuristic mapper
  HeuristicMapper mapper(qc, arch);

  // configure to not include measurements after mapping
  auto config = Configuration{};
  config.addMeasurementsToMappedCircuit = false;

  // perform the mapping
  mapper.map(config);

  // get the resulting circuit
  auto qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper.dumpResult(qasm, qc::Format::OpenQASM3);
  qcMapped.import(qasm, qc::Format::OpenQASM3);

  // check no measurements were added
  EXPECT_EQ(qcMapped.getNops(), 3U);
  EXPECT_NE(qcMapped.back()->getType(), qc::Measure);
}

TEST(Functionality, InvalidCircuits) {
  Configuration config{};
  config.method = Method::Heuristic;

  // architecture not connected
  qc::QuantumComputation qc{2U};
  qc.cx({0}, 1);
  Architecture arch{2U, {}};
  HeuristicMapper mapper(qc, arch);
  EXPECT_THROW(mapper.map(config), QMAPException);

  // gate with >1 control
  qc::QuantumComputation qc3{3U};
  qc3.mcx({0, 1}, 2);
  Architecture arch2{3U, {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 0}, {0, 2}}};
  HeuristicMapper mapper3(qc3, arch2);
  EXPECT_THROW(mapper3.map(config), QMAPException);
}

TEST(Functionality, DataLoggerAfterClose) {
  const std::string dataLoggingPath = "test_log/datalogger_after_close/";
  qc::QuantumComputation qc{3};
  qc.x(0);
  Architecture arch{3, {}};
  auto dataLogger = std::make_unique<DataLogger>(dataLoggingPath, arch, qc);
  const qc::CompoundOperation compOp{};
  Exchange teleport(0, 2, 1, qc::OpType::Teleportation);

  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {{teleport}}, 0);
  dataLogger->logSearchNode(1, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->splitLayer();
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->close();
  dataLogger->clearLog();

  dataLogger->logArchitecture();
  dataLogger->logInputCircuit(qc);
  dataLogger->logOutputCircuit(qc);
  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->splitLayer();
  MappingResults result;
  dataLogger->logMappingResult(result);

  // count files and subdirectories in data logging path
  std::size_t fileCount = 0;
  for ([[maybe_unused]] const auto& _ :
       std::filesystem::directory_iterator(dataLoggingPath)) {
    ++fileCount;
  }
  EXPECT_EQ(fileCount, 0);
}

TEST(Functionality, DataLogger) {
  // setting up example architecture and circuit
  Architecture architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {1, 3}, {3, 1}};
  architecture.loadCouplingMap(4, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.1);
  props.setSingleQubitErrorRate(1, "x", 0.2);
  props.setSingleQubitErrorRate(2, "x", 0.3);
  props.setSingleQubitErrorRate(3, "x", 0.4);

  props.setTwoQubitErrorRate(0, 1, 0.5);
  props.setTwoQubitErrorRate(1, 0, 0.5);
  props.setTwoQubitErrorRate(1, 2, 0.6);
  props.setTwoQubitErrorRate(2, 1, 0.6);
  props.setTwoQubitErrorRate(1, 3, 0.7);
  props.setTwoQubitErrorRate(3, 1, 0.7);

  architecture.loadProperties(props);
  architecture.setName("test_architecture");

  qc::QuantumComputation qc{4, 4};
  qc.cx(1, 0);
  qc.cx(3, 2);
  qc.setName("test_circ");

  Configuration settings{};
  settings.heuristic = Heuristic::GateCountMaxDistance;
  settings.layering = Layering::IndividualGates;
  settings.initialLayout = InitialLayout::Identity;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.lookaheadHeuristic = LookaheadHeuristic::GateCountMaxDistance;
  settings.nrLookaheads = 1;
  settings.firstLookaheadFactor = 0.5;
  settings.lookaheadFactor = 0.9;
  settings.debug = true;
  settings.useTeleportation = false;
  // setting data logging path to enable data logging
  settings.dataLoggingPath = "test_log/datalogger";

  // remove directory at data logging path if it already exists
  if (std::filesystem::exists(settings.dataLoggingPath)) {
    std::filesystem::remove_all(settings.dataLoggingPath);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);
  mapper->map(settings);
  mapper->printResult(std::cout);
  MappingResults& results = mapper->getResults();

  // comparing logged architecture information with original architecture object
  auto archFile =
      std::ifstream(settings.dataLoggingPath + "/architecture.json");
  if (!archFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/architecture.json";
  }
  const auto archJson = nlohmann::basic_json<>::parse(archFile);
  EXPECT_EQ(archJson["name"], architecture.getName());
  EXPECT_EQ(archJson["nqubits"], architecture.getNqubits());
  EXPECT_EQ(archJson["distances"], architecture.getDistanceTable());
  EXPECT_EQ(archJson["coupling_map"], architecture.getCouplingMap());
  const auto& fidelityJson = archJson["fidelity"];
  EXPECT_EQ(fidelityJson["fidelity_distances"],
            architecture.getFidelityDistanceTables());
  EXPECT_EQ(fidelityJson["single_qubit_fidelities"],
            architecture.getSingleQubitFidelities());
  EXPECT_EQ(fidelityJson["two_qubit_fidelities"],
            architecture.getFidelityTable());
  // json does not support inf values, instead nlohmann::json replaces inf
  // with null
  const auto& singleQubitFidelityCosts =
      architecture.getSingleQubitFidelityCosts();
  for (std::size_t i = 0; i < singleQubitFidelityCosts.size(); ++i) {
    EXPECT_EQ(fidelityJson["single_qubit_fidelity_costs"][i],
              singleQubitFidelityCosts[i]);
  }
  const auto& twoQubitFidelityCosts = architecture.getTwoQubitFidelityCosts();
  for (std::size_t i = 0; i < twoQubitFidelityCosts.size(); ++i) {
    for (std::size_t j = 0; j < twoQubitFidelityCosts[i].size(); ++j) {
      EXPECT_EQ(fidelityJson["two_qubit_fidelity_costs"][i][j],
                twoQubitFidelityCosts[i][j]);
    }
  }
  const auto& swapFidelityCosts = architecture.getSwapFidelityCosts();
  for (std::size_t i = 0; i < swapFidelityCosts.size(); ++i) {
    for (std::size_t j = 0; j < swapFidelityCosts[i].size(); ++j) {
      EXPECT_EQ(fidelityJson["swap_fidelity_costs"][i][j],
                swapFidelityCosts[i][j]);
    }
  }

  // comparing logged mapping result with mapping result object
  auto resultFile =
      std::ifstream(settings.dataLoggingPath + "/mapping_result.json");
  if (!resultFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/mapping_result.json";
  }
  const auto resultJson = nlohmann::basic_json<>::parse(resultFile);
  const auto& configJson = resultJson["config"];
  EXPECT_EQ(configJson["add_measurements_to_mapped_circuit"],
            settings.addMeasurementsToMappedCircuit);
  EXPECT_EQ(configJson["debug"], settings.debug);
  EXPECT_EQ(configJson["verbose"], settings.verbose);
  EXPECT_EQ(configJson["layering_strategy"], toString(settings.layering));
  EXPECT_EQ(configJson["method"], toString(settings.method));
  EXPECT_EQ(configJson["post_mapping_optimizations"],
            settings.postMappingOptimizations);
  EXPECT_EQ(configJson["pre_mapping_optimizations"],
            settings.preMappingOptimizations);
  const auto& heuristicSettingsJson = configJson["settings"];
  EXPECT_EQ(heuristicSettingsJson["heuristic"], toString(settings.heuristic));
  EXPECT_EQ(heuristicSettingsJson["initial_layout"],
            toString(settings.initialLayout));
  if (settings.lookaheadHeuristic != LookaheadHeuristic::None) {
    const auto& lookaheadJson = heuristicSettingsJson["lookahead"];
    EXPECT_EQ(lookaheadJson["heuristic"],
              toString(settings.lookaheadHeuristic));
    EXPECT_EQ(lookaheadJson["factor"], settings.lookaheadFactor);
    EXPECT_EQ(lookaheadJson["first_factor"], settings.firstLookaheadFactor);
    EXPECT_EQ(lookaheadJson["lookaheads"], settings.nrLookaheads);
  }
  if (settings.useTeleportation) {
    const auto& teleportationJson = heuristicSettingsJson["teleportation"];
    EXPECT_EQ(teleportationJson["qubits"], settings.teleportationQubits);
    EXPECT_EQ(teleportationJson["seed"], settings.teleportationSeed);
    EXPECT_EQ(teleportationJson["fake"], settings.teleportationFake);
  }

  const auto& inCircJson = resultJson["circuit"];
  EXPECT_EQ(inCircJson["cnots"], results.input.cnots);
  EXPECT_EQ(inCircJson["gates"], results.input.gates);
  EXPECT_EQ(inCircJson["name"], results.input.name);
  EXPECT_EQ(inCircJson["qubits"], results.input.qubits);
  EXPECT_EQ(inCircJson["single_qubit_gates"], results.input.singleQubitGates);

  const auto& outCircJson = resultJson["mapped_circuit"];
  EXPECT_EQ(outCircJson["cnots"], results.output.cnots);
  EXPECT_EQ(outCircJson["gates"], results.output.gates);
  EXPECT_EQ(outCircJson["name"], results.output.name);
  EXPECT_EQ(outCircJson["qubits"], results.output.qubits);
  EXPECT_EQ(outCircJson["single_qubit_gates"], results.output.singleQubitGates);

  const auto& statJson = resultJson["statistics"];
  EXPECT_EQ(statJson["additional_gates"],
            (static_cast<std::make_signed_t<decltype(results.output.gates)>>(
                 results.output.gates) -
             static_cast<std::make_signed_t<decltype(results.input.gates)>>(
                 results.input.gates)));
  EXPECT_EQ(statJson["layers"], results.input.layers);
  EXPECT_EQ(statJson["mapping_time"], results.time);
  EXPECT_EQ(statJson["swaps"], results.output.swaps);
  EXPECT_EQ(statJson["teleportations"], results.output.teleportations);
  EXPECT_EQ(statJson["total_fidelity"], results.output.totalFidelity);
  EXPECT_EQ(statJson["total_log_fidelity"], results.output.totalLogFidelity);

  const auto& benchmarkJson = statJson["benchmark"];
  EXPECT_EQ(benchmarkJson["average_branching_factor"],
            results.heuristicBenchmark.averageBranchingFactor);
  EXPECT_EQ(benchmarkJson["effective_branching_factor"],
            results.heuristicBenchmark.effectiveBranchingFactor);
  EXPECT_EQ(benchmarkJson["expanded_nodes"],
            results.heuristicBenchmark.expandedNodes);
  EXPECT_EQ(benchmarkJson["generated_nodes"],
            results.heuristicBenchmark.generatedNodes);
  EXPECT_EQ(benchmarkJson["seconds_per_node"],
            results.heuristicBenchmark.secondsPerNode);
  const auto& benchmarkLayersJson = benchmarkJson["layers"];
  EXPECT_EQ(benchmarkLayersJson.size(), results.layerHeuristicBenchmark.size());
  for (std::size_t i = 0; i < results.layerHeuristicBenchmark.size(); ++i) {
    EXPECT_EQ(benchmarkLayersJson[i]["average_branching_factor"],
              results.layerHeuristicBenchmark.at(i).averageBranchingFactor);
    EXPECT_EQ(benchmarkLayersJson[i]["effective_branching_factor"],
              results.layerHeuristicBenchmark.at(i).effectiveBranchingFactor);
    EXPECT_EQ(benchmarkLayersJson[i]["expanded_nodes"],
              results.layerHeuristicBenchmark.at(i).expandedNodes);
    EXPECT_EQ(benchmarkLayersJson[i]["generated_nodes"],
              results.layerHeuristicBenchmark.at(i).generatedNodes);
    EXPECT_EQ(benchmarkLayersJson[i]["solution_depth"],
              results.layerHeuristicBenchmark.at(i).solutionDepth);
    EXPECT_EQ(benchmarkLayersJson[i]["seconds_per_node"],
              results.layerHeuristicBenchmark.at(i).secondsPerNode);
  }

  // comparing logged input and output circuits with input circuit object and
  // mapped circuit object
  auto inputQasmFile = std::ifstream(settings.dataLoggingPath + "/input.qasm");
  if (!inputQasmFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/input.qasm";
  }
  std::stringstream inputFileBuffer;
  inputFileBuffer << inputQasmFile.rdbuf();
  std::stringstream inputQasmBuffer;
  qc.dumpOpenQASM3(inputQasmBuffer);
  EXPECT_EQ(inputFileBuffer.str(), inputQasmBuffer.str());

  auto outputQasmFile =
      std::ifstream(settings.dataLoggingPath + "/output.qasm");
  if (!outputQasmFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/output.qasm";
  }
  std::stringstream outputFileBuffer;
  outputFileBuffer << outputQasmFile.rdbuf();
  std::stringstream outputQasmBuffer;
  mapper->dumpResult(outputQasmBuffer, qc::Format::OpenQASM3);
  EXPECT_EQ(outputFileBuffer.str(), outputQasmBuffer.str());

  // checking logged search graph info against known values (correct qubit
  // number, valid layouts, correct data types in all csv fields, etc.)
  for (std::size_t i = 0; i < results.input.layers; ++i) {
    auto layerFile = std::ifstream(settings.dataLoggingPath + "/layer_" +
                                   std::to_string(i) + ".json");
    if (!layerFile.is_open()) {
      FAIL() << "Could not open file " << settings.dataLoggingPath << "/layer_"
             << i << ".json";
    }
    const auto layerJson = nlohmann::basic_json<>::parse(layerFile);
    const std::size_t finalNodeId = layerJson["final_node_id"];
    EXPECT_EQ(layerJson["initial_layout"].size(), architecture.getNqubits());
    EXPECT_EQ(layerJson["single_qubit_multiplicity"].size(),
              architecture.getNqubits());

    std::vector<HeuristicMapper::Node> nodes{
        results.layerHeuristicBenchmark.at(i).generatedNodes,
        HeuristicMapper::Node{architecture.getNqubits(), 0}};
    parseNodesFromDatalog(settings.dataLoggingPath, i, nodes);

    if (finalNodeId >= nodes.size() ||
        nodes.at(finalNodeId).id != finalNodeId) {
      FAIL() << "Final solution node " << finalNodeId
             << " not found in nodes of layer " << i;
    }
    auto& finalSolutionNode = nodes.at(finalNodeId);
    EXPECT_EQ(layerJson["final_cost_fixed"], finalSolutionNode.costFixed);
    EXPECT_EQ(layerJson["final_cost_heur"], finalSolutionNode.costHeur);
    EXPECT_EQ(layerJson["final_lookahead_penalty"],
              finalSolutionNode.lookaheadPenalty);
    EXPECT_EQ(layerJson["final_search_depth"], finalSolutionNode.depth);
    std::vector<std::int16_t> layout{};
    for (std::size_t j = 0; j < architecture.getNqubits(); ++j) {
      layout.emplace_back(finalSolutionNode.qubits.at(j));
    }
    std::vector<std::pair<std::uint16_t, std::uint16_t>> swaps{};
    swaps.reserve(finalSolutionNode.swaps.size());
    for (auto& swap : finalSolutionNode.swaps) {
      swaps.emplace_back(swap.first, swap.second);
    }
    EXPECT_EQ(layerJson["final_layout"], layout);
    EXPECT_EQ(layerJson["final_swaps"], swaps);
    EXPECT_EQ(finalSolutionNode.validMapping, true);

    for (std::size_t j = 0; j < nodes.size(); ++j) {
      auto& node = nodes.at(j);
      if (j != node.id) {
        continue;
      }

      for (std::size_t k = 0; k < architecture.getNqubits(); ++k) {
        EXPECT_GE(node.qubits.at(k), -1);
        EXPECT_LT(node.qubits.at(k), qc.getNqubits());
        EXPECT_GE(node.locations.at(k), -1);
        EXPECT_LT(node.locations.at(k), architecture.getNqubits());

        if (node.qubits.at(k) >= 0) {
          EXPECT_EQ(
              node.locations.at(static_cast<std::size_t>(node.qubits.at(k))),
              static_cast<std::int16_t>(k));
        }
        if (node.locations.at(k) >= 0) {
          EXPECT_EQ(
              node.qubits.at(static_cast<std::size_t>(node.locations.at(k))),
              static_cast<std::int16_t>(k));
        }
      }
    }
  }

  // checking if files for non-existing layers are not created
  auto afterLastLayerFile =
      std::ifstream(settings.dataLoggingPath + "/layer_" +
                    std::to_string(results.input.layers) + ".json");
  if (afterLastLayerFile.is_open()) {
    FAIL() << "File " << settings.dataLoggingPath << "/layer_"
           << results.input.layers
           << ".json should not exist, as there are not that many layers";
  }
  auto afterLastLayerNodesFile =
      std::ifstream(settings.dataLoggingPath + "/nodes_layer_" +
                    std::to_string(results.input.layers) + ".csv");
  if (afterLastLayerNodesFile.is_open()) {
    FAIL() << "File " << settings.dataLoggingPath << "/nodes_layer_"
           << results.input.layers
           << ".csv should not exist, as there are not that many layers";
  }
}

TEST(Functionality, terminationStrategyFromString) {
  const std::vector<std::pair<std::string, EarlyTermination>>
      terminationStrategies = {
          {"none", EarlyTermination::None},
          {"expanded_nodes", EarlyTermination::ExpandedNodes},
          {"expanded_nodes_after_first_solution",
           EarlyTermination::ExpandedNodesAfterFirstSolution},
          {"expanded_nodes_after_current_optimal_solution",
           EarlyTermination::ExpandedNodesAfterCurrentOptimalSolution},
          {"solution_nodes", EarlyTermination::SolutionNodes},
          {"solution_nodes_after_current_optimal_solution",
           EarlyTermination::SolutionNodesAfterCurrentOptimalSolution}};

  for (const auto& [str, termination] : terminationStrategies) {
    EXPECT_EQ(earlyTerminationFromString(str), termination);
  }
  EXPECT_THROW(earlyTerminationFromString("invalid"), std::invalid_argument);
}

TEST(Functionality, earlyTermination) {
  qc::QuantumComputation qc{7, 7};
  qc.x(0);
  qc.cx(1, 2);
  for (std::size_t i = 0; i < 7; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                          {3, 4}, {4, 3}, {4, 5}, {5, 4}, {5, 6}, {6, 5}};
  Architecture::Properties props{};
  props.setSingleQubitErrorRate(0, "x", 0.9);
  props.setSingleQubitErrorRate(1, "x", 0.5);
  props.setSingleQubitErrorRate(2, "x", 0.5);
  props.setSingleQubitErrorRate(3, "x", 0.5);
  props.setSingleQubitErrorRate(4, "x", 0.5);
  props.setSingleQubitErrorRate(5, "x", 0.5);
  props.setSingleQubitErrorRate(6, "x", 0.1);
  for (auto edge : cm) {
    props.setTwoQubitErrorRate(edge.first, edge.second, 0.01, "cx");
  }
  Architecture arch{7, cm, props};

  Configuration config{};
  config.method = Method::Heuristic;
  config.layering = Layering::DisjointQubits;
  config.initialLayout = InitialLayout::Identity;
  config.automaticLayerSplits = false;
  config.iterativeBidirectionalRouting = false;
  config.debug = true;
  config.heuristic = Heuristic::FidelityBestLocation;
  config.lookaheadHeuristic = LookaheadHeuristic::None;
  config.earlyTerminationLimit = 4;

  auto mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination = EarlyTermination::None;
  mapper->map(config);
  auto results = mapper->getResults();
  EXPECT_FALSE(results.layerHeuristicBenchmark[0].earlyTermination);

  mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination = EarlyTermination::ExpandedNodesAfterFirstSolution;
  mapper->map(config);
  results = mapper->getResults();
  EXPECT_TRUE(results.layerHeuristicBenchmark[0].earlyTermination);
  EXPECT_EQ(results.layerHeuristicBenchmark[0].expandedNodesAfterFirstSolution,
            4);

  mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination =
      EarlyTermination::ExpandedNodesAfterCurrentOptimalSolution;
  mapper->map(config);
  results = mapper->getResults();
  EXPECT_TRUE(results.layerHeuristicBenchmark[0].earlyTermination);
  EXPECT_EQ(
      results.layerHeuristicBenchmark[0].expandedNodesAfterOptimalSolution, 4);

  mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination = EarlyTermination::ExpandedNodes;
  mapper->map(config);
  results = mapper->getResults();
  EXPECT_TRUE(results.layerHeuristicBenchmark[0].earlyTermination);
  EXPECT_EQ(results.layerHeuristicBenchmark[0].expandedNodes, 4);

  mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination = EarlyTermination::SolutionNodes;
  mapper->map(config);
  results = mapper->getResults();
  EXPECT_TRUE(results.layerHeuristicBenchmark[0].earlyTermination);
  EXPECT_EQ(results.layerHeuristicBenchmark[0].solutionNodes, 4);

  mapper = std::make_unique<HeuristicMapper>(qc, arch);
  config.earlyTermination =
      EarlyTermination::SolutionNodesAfterCurrentOptimalSolution;
  mapper->map(config);
  results = mapper->getResults();
  EXPECT_TRUE(results.layerHeuristicBenchmark[0].earlyTermination);
  EXPECT_EQ(
      results.layerHeuristicBenchmark[0].solutionNodesAfterOptimalSolution, 4);
}

TEST(Functionality, InitialLayoutDump) {
  // queko's BNTF/16QBT_05CYC_TFL_9.qasm
  qc::QuantumComputation qc{16U};
  qc.x(15);
  qc.x(4);
  qc.cx(10, 8);
  qc.cx(3, 5);
  qc.cx(7, 1);
  qc.cx(6, 12);
  qc.cx(0, 13);
  qc.x(10);
  qc.x(9);
  qc.x(1);
  qc.x(2);
  qc.x(12);
  qc.cx(3, 11);
  qc.cx(0, 13);
  qc.x(9);
  qc.x(6);
  qc.x(2);
  qc.x(0);
  qc.x(12);
  qc.x(7);
  qc.cx(10, 8);
  qc.cx(13, 5);
  qc.x(10);
  qc.x(15);
  qc.x(9);
  qc.x(1);
  qc.x(13);
  qc.x(5);
  qc.cx(6, 12);
  qc.cx(7, 11);
  qc.cx(14, 4);
  qc.x(1);
  qc.x(14);
  qc.x(4);
  qc.cx(12, 10);
  qc.cx(3, 5);
  qc.cx(8, 7);

  // subgraph of IBM's Brisbane Backend
  Architecture arch{27U,
                    {{1, 0},   {2, 1},   {3, 2},   {4, 3},   {4, 5},   {4, 15},
                     {6, 5},   {6, 7},   {7, 8},   {8, 9},   {10, 9},  {10, 11},
                     {11, 12}, {12, 17}, {13, 12}, {14, 0},  {14, 18}, {15, 22},
                     {16, 8},  {16, 26}, {18, 19}, {20, 19}, {21, 20}, {21, 22},
                     {22, 23}, {24, 23}, {25, 24}, {26, 25}}};

  Configuration config{};
  config.method = Method::Heuristic;
  config.layering = Layering::Disjoint2qBlocks;

  HeuristicMapper mapper(qc, arch);
  mapper.map(config);

  std::stringstream qasmStream{};
  mapper.dumpResult(qasmStream, qc::Format::OpenQASM3);
  const std::string qasm = qasmStream.str();

  qasmStream = std::stringstream(qasm);
  auto qcMapped = qc::QuantumComputation();
  qcMapped.import(qasmStream, qc::Format::OpenQASM3);

  qasmStream = std::stringstream(qasm);
  std::string line;
  bool foundPermutation = false;
  while (std::getline(qasmStream, line)) {
    if (line.rfind("// i ", 0) == 0) {
      std::stringstream lineStream(line.substr(5));
      std::string entry;
      std::vector<std::uint32_t> qubits{};
      while (std::getline(lineStream, entry, ' ')) {
        EXPECT_NO_THROW(
            qubits.emplace_back(static_cast<std::uint32_t>(std::stoul(entry))))
            << "invalid qubit id " << entry;
      }
      const std::set<std::uint32_t> qubitSet(qubits.begin(), qubits.end());
      EXPECT_EQ(qubitSet.size(), qubits.size());
      for (std::uint32_t i = 0; i < qcMapped.getNqubits(); ++i) {
        EXPECT_TRUE(qubitSet.count(i) > 0)
            << "qubit " << std::to_string(i) << " not found in layout";
      }
      foundPermutation = true;
      break;
    }
  }
  EXPECT_TRUE(foundPermutation) << "no initial layout found in mapped circuit";
}

class LayeringTest : public testing::Test {
protected:
  qc::QuantumComputation qc;
  Architecture arch;
  std::unique_ptr<HeuristicMapper> mapper;
  Configuration settings{};

  void SetUp() override {
    qc = qc::QuantumComputation{4, 4};
    qc.x(0);
    qc.x(1);
    qc.cx(0, 1);
    qc.cx(2, 3);
    qc.cx(1, 2);
    qc.x(3);
    qc.barrier({0, 1, 2});
    for (size_t i = 0; i < 3; ++i) {
      qc.measure(static_cast<qc::Qubit>(i), i);
    }

    arch = Architecture{4, {{0, 1}, {1, 2}, {2, 3}}};

    settings.initialLayout = InitialLayout::Dynamic;
    settings.preMappingOptimizations = false;
    settings.postMappingOptimizations = false;
    settings.addMeasurementsToMappedCircuit = true;
    settings.addBarriersBetweenLayers = true;
    settings.automaticLayerSplits = false;

    mapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

TEST_F(LayeringTest, Disjoint2qBlocks) {
  settings.layering = Layering::Disjoint2qBlocks;
  mapper->map(settings);
  auto result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 2);
  // get mapped circuit
  auto qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper->dumpResult(qasm, qc::Format::OpenQASM3);
  qcMapped.import(qasm, qc::Format::OpenQASM3);
  // check barrier count
  std::size_t barriers = 0;
  for (const auto& op : qcMapped) {
    if (op->getType() == qc::Barrier) {
      ++barriers;
    }
  }
  EXPECT_EQ(barriers, result.input.layers);
}

TEST_F(LayeringTest, DisjointQubits) {
  settings.layering = Layering::DisjointQubits;
  mapper->map(settings);
  auto result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 3);
  // get mapped circuit
  auto qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper->dumpResult(qasm, qc::Format::OpenQASM3);
  qcMapped.import(qasm, qc::Format::OpenQASM3);
  // check barrier count
  std::size_t barriers = 0;
  for (const auto& op : qcMapped) {
    if (op->getType() == qc::Barrier) {
      ++barriers;
    }
  }
  EXPECT_EQ(barriers, result.input.layers);
}

TEST_F(LayeringTest, IndividualGates) {
  settings.layering = Layering::IndividualGates;
  mapper->map(settings);
  auto result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 6);
  // get mapped circuit
  auto qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper->dumpResult(qasm, qc::Format::OpenQASM3);
  qcMapped.import(qasm, qc::Format::OpenQASM3);
  // check barrier count
  std::size_t barriers = 0;
  for (const auto& op : qcMapped) {
    if (op->getType() == qc::Barrier) {
      ++barriers;
    }
  }
  EXPECT_EQ(barriers, result.input.layers);
}

class HeuristicTest5Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";
  std::string testCalibrationDir = "../../../extern/calibration/";

  qc::QuantumComputation qc;
  Architecture ibmqYorktown;
  Architecture ibmqLondon;
  std::unique_ptr<HeuristicMapper> ibmqYorktownMapper;
  std::unique_ptr<HeuristicMapper> ibmqLondonMapper;
  Configuration settings{};

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmqYorktown.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
    ibmqYorktownMapper = std::make_unique<HeuristicMapper>(qc, ibmqYorktown);
    ibmqLondonMapper = std::make_unique<HeuristicMapper>(qc, ibmqLondon);
    settings.verbose = true;
    settings.debug = true;

    settings.iterativeBidirectionalRouting = true;
    settings.iterativeBidirectionalRoutingPasses = 3;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest5Q,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest5Q, Identity) {
  settings.initialLayout = InitialLayout::Identity;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
  settings.initialLayout = InitialLayout::Static;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
  settings.initialLayout = InitialLayout::Dynamic;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";

  qc::QuantumComputation qc;
  Architecture ibmQX5;
  std::unique_ptr<HeuristicMapper> ibmQX5Mapper;
  Configuration settings{};

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmQX5.loadCouplingMap(AvailableArchitecture::IbmQx5);
    ibmQX5Mapper = std::make_unique<HeuristicMapper>(qc, ibmQX5);
    settings.debug = true;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest16Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16"),
    [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest16Q, Dynamic) {
  settings.initialLayout = InitialLayout::Dynamic;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint) {
  settings.layering = Layering::DisjointQubits;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint2qBlocks) {
  settings.layering = Layering::Disjoint2qBlocks;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";

  qc::QuantumComputation qc;
  Architecture arch;
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest20Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16",
                    "z4_268"),
    [](const testing::TestParamInfo<HeuristicTest20Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest20Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  settings.debug = true;
  tokyoMapper->map(settings);
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20QTeleport
    : public testing::TestWithParam<std::tuple<std::uint64_t, std::string>> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";

  qc::QuantumComputation qc;
  Architecture arch;
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + std::get<1>(GetParam()) + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    HeuristicTeleport, HeuristicTest20QTeleport,
    testing::Combine(testing::Values(0, 1, 2, 3, 1337, 1338, 3147),
                     testing::Values("ising_model_10", "rd73_140", "cnt3-5_179",
                                     "qft_16", "z4_268")),
    [](const testing::TestParamInfo<HeuristicTest20QTeleport::ParamType>& inf) {
      std::string name = std::get<1>(inf.param);
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name << "_seed" << std::get<0>(inf.param);
      return ss.str();
    });

TEST_P(HeuristicTest20QTeleport, Teleportation) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  settings.debug = true;
  settings.teleportationQubits = std::min(
      (arch.getNqubits() - qc.getNqubits()) & ~1U, static_cast<std::size_t>(8));
  settings.teleportationSeed = std::get<0>(GetParam());
  tokyoMapper->map(settings);
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTestFidelity : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir = "../../../examples/";
  std::string testArchitectureDir = "../../../extern/architectures/";
  std::string testCalibrationDir = "../../../extern/calibration/";

  qc::QuantumComputation qc;
  Architecture arch;
  Architecture nonFidelityArch;
  std::unique_ptr<HeuristicMapper> mapper;
  std::unique_ptr<HeuristicMapper> nonFidelityMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    arch.loadProperties(testCalibrationDir + "ibmq_london.csv");
    mapper = std::make_unique<HeuristicMapper>(qc, arch);
    nonFidelityArch.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    nonFidelityMapper = std::make_unique<HeuristicMapper>(qc, nonFidelityArch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTestFidelity,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTestFidelity::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTestFidelity, Identity) {
  Configuration settings{};
  settings.layering = Layering::DisjointQubits;
  settings.initialLayout = InitialLayout::Identity;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  mapper->map(settings);
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, Static) {
  Configuration settings{};
  settings.layering = Layering::DisjointQubits;
  settings.initialLayout = InitialLayout::Static;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  mapper->map(settings);
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, Dynamic) {
  Configuration settings{};
  settings.layering = Layering::DisjointQubits;
  settings.initialLayout = InitialLayout::Dynamic;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  mapper->map(settings);
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, NoFidelity) {
  Configuration settings{};
  settings.layering = Layering::DisjointQubits;
  settings.initialLayout = InitialLayout::Static;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  EXPECT_THROW(nonFidelityMapper->map(settings), QMAPException);
}

TEST(HeuristicTestFidelity, RemapSingleQubit) {
  Architecture architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
      {3, 2}, {3, 4}, {4, 3}, {4, 5}, {5, 4},
  };
  architecture.loadCouplingMap(6, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e1 = 0.1;
  const double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e4);
  props.setSingleQubitErrorRate(4, "x", e4);
  props.setSingleQubitErrorRate(5, "x", e1);

  props.setTwoQubitErrorRate(0, 1, e1);
  props.setTwoQubitErrorRate(1, 0, e1);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e5);
  props.setTwoQubitErrorRate(3, 2, e5);
  props.setTwoQubitErrorRate(3, 4, e0);
  props.setTwoQubitErrorRate(4, 3, e0);
  props.setTwoQubitErrorRate(4, 5, e0);
  props.setTwoQubitErrorRate(5, 4, e0);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{6, 6};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(2, 0);
    qc.x(3);
  }

  for (size_t i = 0; i < 6; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.layering = Layering::Disjoint2qBlocks;
  settings.initialLayout = InitialLayout::Identity;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.verbose = true;
  settings.swapOnFirstLayer = true;
  mapper->map(settings);
  mapper->printResult(std::cout);

  // 0 --e1-- 1 --e3-- 2 --e5-- 3 --e0-- 4 --e0-- 5
  // e5       e5       e5       e4       e4       e1

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(1,2)
  SWAP(3,4)
  SWAP(4,5)
  CX(0,1) [x5]
  X(5) [x5]
  */
  EXPECT_EQ(result.output.swaps, 3);

  const double c3 = -std::log2(1 - e3);
  const double c1 = -std::log2(1 - e1);
  const double c0 = -std::log2(1 - e0);

  const double expectedFidelity = 3 * c3 + 3 * c0 + 3 * c0 + // SWAPs
                                  5 * c1 +                   // Xs
                                  5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, QubitRideAlong) {
  Architecture architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                          {1, 4}, {4, 1}, {2, 5}, {5, 2}, {5, 6}, {6, 5}};
  architecture.loadCouplingMap(7, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e1 = 0.1;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e5);
  props.setSingleQubitErrorRate(5, "x", e5);
  props.setSingleQubitErrorRate(6, "x", e5);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);
  props.setTwoQubitErrorRate(1, 4, e1);
  props.setTwoQubitErrorRate(4, 1, e1);
  props.setTwoQubitErrorRate(2, 5, e3);
  props.setTwoQubitErrorRate(5, 2, e3);
  props.setTwoQubitErrorRate(5, 6, e4);
  props.setTwoQubitErrorRate(6, 5, e4);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{7, 7};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(3, 0);
    qc.cx(6, 4);
  }

  for (size_t i = 0; i < 7; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.layering = Layering::Disjoint2qBlocks;
  settings.initialLayout = InitialLayout::Identity;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.verbose = true;
  settings.swapOnFirstLayer = true;
  mapper->map(settings);
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(5,6)
  SWAP(2,5)
  SWAP(0,1)
  SWAP(1,2)
  CX(2,3) [x5]
  CX(4,1) [x5]
  */
  EXPECT_EQ(result.output.swaps, 4);

  const double c4 = -std::log2(1 - e4);
  const double c3 = -std::log2(1 - e3);
  const double c1 = -std::log2(1 - e1);

  const double expectedFidelity = 3 * c4 + 3 * c3 + 3 * c4 + 3 * c3 + // SWAPs
                                  5 * c1 + 5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, SingleQubitsCompete) {
  Architecture architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}};
  architecture.loadCouplingMap(3, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.8);
  props.setSingleQubitErrorRate(1, "x", 0.1);
  props.setSingleQubitErrorRate(2, "x", 0.8);

  props.setTwoQubitErrorRate(0, 1, 0.1);
  props.setTwoQubitErrorRate(1, 0, 0.1);
  props.setTwoQubitErrorRate(1, 2, 0.1);
  props.setTwoQubitErrorRate(2, 1, 0.1);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{3, 3};
  qc.x(0);
  qc.x(2);

  for (size_t i = 0; i < 3; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.layering = Layering::Disjoint2qBlocks;
  settings.initialLayout = InitialLayout::Identity;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.verbose = true;
  settings.swapOnFirstLayer = true;
  mapper->map(settings);
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  EXPECT_EQ(result.output.swaps, 1);

  const double expectedFidelity = -3 * std::log2(1 - 0.1) // SWAPs
                                  - std::log2(1 - 0.8) -
                                  std::log2(1 - 0.1); // Xs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, LayerSplitting) {
  Architecture architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2},  {2, 1},  {2, 3},   {3, 2},

      {0, 4}, {4, 0}, {1, 5},  {5, 1},  {2, 6},   {6, 2},  {3, 7},  {7, 3},

      {4, 5}, {5, 4}, {5, 6},  {6, 5},  {6, 7},   {7, 6},

      {4, 8}, {8, 4}, {5, 9},  {9, 5},  {6, 10},  {10, 6}, {7, 11}, {11, 7},

      {8, 9}, {9, 8}, {9, 10}, {10, 9}, {10, 11}, {11, 10}};
  architecture.loadCouplingMap(12, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e2 = 0.4;
  const double e1 = 0.1;
  const double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e3);
  props.setSingleQubitErrorRate(5, "x", e3);
  props.setSingleQubitErrorRate(6, "x", e5);
  props.setSingleQubitErrorRate(7, "x", e3);
  props.setSingleQubitErrorRate(8, "x", e2);
  props.setSingleQubitErrorRate(9, "x", e1);
  props.setSingleQubitErrorRate(10, "x", e5);
  props.setSingleQubitErrorRate(11, "x", e2);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e4);
  props.setTwoQubitErrorRate(2, 1, e4);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);

  props.setTwoQubitErrorRate(0, 4, e5);
  props.setTwoQubitErrorRate(4, 0, e5);
  props.setTwoQubitErrorRate(1, 5, e5);
  props.setTwoQubitErrorRate(5, 1, e5);
  props.setTwoQubitErrorRate(2, 6, e4);
  props.setTwoQubitErrorRate(6, 2, e4);
  props.setTwoQubitErrorRate(3, 7, e5);
  props.setTwoQubitErrorRate(7, 3, e5);

  props.setTwoQubitErrorRate(4, 5, e3);
  props.setTwoQubitErrorRate(5, 4, e3);
  props.setTwoQubitErrorRate(5, 6, e5);
  props.setTwoQubitErrorRate(6, 5, e5);
  props.setTwoQubitErrorRate(6, 7, e5);
  props.setTwoQubitErrorRate(7, 6, e5);

  props.setTwoQubitErrorRate(4, 8, e0);
  props.setTwoQubitErrorRate(8, 4, e0);
  props.setTwoQubitErrorRate(5, 9, e3);
  props.setTwoQubitErrorRate(9, 5, e3);
  props.setTwoQubitErrorRate(6, 10, e1);
  props.setTwoQubitErrorRate(10, 6, e1);
  props.setTwoQubitErrorRate(7, 11, e3);
  props.setTwoQubitErrorRate(11, 7, e3);

  props.setTwoQubitErrorRate(8, 9, e5);
  props.setTwoQubitErrorRate(9, 8, e5);
  props.setTwoQubitErrorRate(9, 10, e4);
  props.setTwoQubitErrorRate(10, 9, e4);
  props.setTwoQubitErrorRate(10, 11, e5);
  props.setTwoQubitErrorRate(11, 10, e5);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{12, 12};

  for (std::size_t i = 0; i < 50; ++i) {
    qc.x(4);
  }
  qc.x(5);
  qc.x(7);
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(3, 0);
    qc.cx(9, 2);
  }

  for (size_t i = 0; i < 12; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.verbose = true;
  settings.layering = Layering::Disjoint2qBlocks;
  settings.initialLayout = InitialLayout::Identity;
  settings.heuristic = Heuristic::FidelityBestLocation;
  settings.lookaheadHeuristic = LookaheadHeuristic::None;
  settings.preMappingOptimizations = false;
  settings.postMappingOptimizations = false;
  settings.swapOnFirstLayer = true;
  settings.automaticLayerSplits = true;
  settings.automaticLayerSplitsNodeLimit =
      1; // force splittings after 1st expanded node until layers are
         // unsplittable
  settings.dataLoggingPath = "test_log/layer_splitting/";
  mapper->map(settings);
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 5); // originally 1 but split into 5 during A*
  /*
  expected output:
  === layout (placing of logical qubits in physical grid):
       0  1  2  3
       4  5  6  7
       8  9 10 11
  SWAP (4,5)
  SWAP (5,9)
  === layout:
       0  1  2  3
       5  9  6  7
       8  4 10 11
  X(9) [x50]   (originally X(4))

  === layout:
       0  1  2  3
       5  9  6  7
       8  4 10 11
  X(7)   (originally X(7))

  SWAP (5,9)
  SWAP (9,10)
  SWAP (2,6)
  === layout:
       0  1  6  3
       5  4  2  7
       8 10  9 11
  CX(6,10) [x5]   (originally CX(2,9))

  SWAP (4,8)
  === layout:
       0  1  6  3
       8  4  2  7
       5 10  9 11
  X(8)   (originally X(5))

  SWAP (0,1)
  SWAP (1,2)
  === layout:
       1  6  0  3
       8  4  2  7
       5 10  9 11
  CX(2,3) [x5]   (originally CX(0,3))
  */
  EXPECT_EQ(result.output.swaps, 8);

  const double c4 = -std::log2(1 - e4);
  const double c3 = -std::log2(1 - e3);
  const double c2 = -std::log2(1 - e2);
  const double c1 = -std::log2(1 - e1);
  const double c0 = -std::log2(1 - e0);

  const double expectedFidelity =
      3 * (c3 + c3 + c3 + c4 + c4 + c0 + c4 + c4) + // SWAPs
      50 * c1 + c3 + c2 +                           // Xs
      5 * c1 + 5 * c1;                              // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);

  // check data log
  const std::array<std::string, 4> layerNodeFilePaths = {
      "nodes_layer_0.presplit-0.csv", "nodes_layer_0.presplit-1.csv",
      "nodes_layer_1.presplit-0.csv", "nodes_layer_3.presplit-0.csv"};
  for (const auto& path : layerNodeFilePaths) {
    auto layerNodeFile = std::ifstream(settings.dataLoggingPath + path);
    if (!layerNodeFile.is_open()) {
      FAIL() << "Could not open file " << settings.dataLoggingPath << path;
    }
    std::string line;
    while (std::getline(layerNodeFile, line)) {
      if (line.empty()) {
        continue;
      }
      std::string col;
      std::stringstream lineStream(line);
      if (!std::getline(lineStream, col, ';')) {
        FAIL() << "Missing value for node id in " << settings.dataLoggingPath
               << path;
      }
      if (std::getline(lineStream, col, ';')) {
        if (std::stoull(col) != 0) {
          // should only contain root node and its direct children
          FAIL() << "Unexpected value for parent node id in "
                 << settings.dataLoggingPath << path;
        }
      } else {
        FAIL() << "Missing value for parent node id in "
               << settings.dataLoggingPath << path;
      }
    }
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_0.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_0.presplit-0.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_0.presplit-1.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_0.presplit-1.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_1.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_1.presplit-0.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_3.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_3.presplit-0.json"
           << " does not exist";
  }
}

TEST(HeuristicDebug, MoreThan128Qubits) {
  const std::size_t numQubits = 129;

  // create a linear NN coupling map
  CouplingMap cm = {};
  for (std::size_t i = 0; i < numQubits - 1; ++i) {
    cm.emplace(i, i + 1);
    cm.emplace(i + 1, i);
  }
  Architecture architecture{numQubits, cm};

  // create a simple circuit to map to the architecture
  auto qc = qc::QuantumComputation(numQubits);
  for (std::size_t i = 0; i < numQubits - 1; ++i) {
    qc.cx(static_cast<qc::Qubit>(i), static_cast<qc::Qubit>(i + 1));
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};

  settings.verbose = true;
  mapper->map(settings);
  mapper->printResult(std::cout);
}
