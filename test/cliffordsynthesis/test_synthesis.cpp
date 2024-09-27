//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Definitions.hpp"
#include "cliffordsynthesis/CliffordSynthesizer.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/TargetMetric.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Control.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>
#include <plog/Severity.h>
#include <sstream>
#include <string>
#include <vector>

using namespace qc::literals;

namespace cs {

struct TestConfiguration {
  // given input (either as tableau or as circuit)
  std::string description;
  std::string initialTableau;
  std::string targetTableau;
  std::string initialCircuit;
  std::string couplingMap;

  // expected output
  std::size_t expectedMinimalGates{};
  std::size_t expectedMinimalDepth{};
  std::size_t expectedMinimalGatesAtMinimalDepth{};
  std::size_t expectedMinimalTwoQubitGates{};
  std::size_t expectedMinimalGatesAtMinimalTwoQubitGates{};
};

// NOLINTNEXTLINE (readability-identifier-naming)
inline void from_json(const nlohmann::json& j, TestConfiguration& test) {
  test.description = j.at("description").get<std::string>();
  if (j.contains("initial_tableau")) {
    test.initialTableau = j.at("initial_tableau").get<std::string>();
  }
  if (j.contains("target_tableau")) {
    test.targetTableau = j.at("target_tableau").get<std::string>();
  }
  if (j.contains("initial_circuit")) {
    test.initialCircuit = j.at("initial_circuit").get<std::string>();
  }
  if (j.contains("coupling_map")) {
    test.couplingMap = j.at("coupling_map").get<std::string>();
  }

  test.expectedMinimalGates = j.at("expected_minimal_gates").get<std::size_t>();
  test.expectedMinimalDepth = j.at("expected_minimal_depth").get<std::size_t>();
  test.expectedMinimalGatesAtMinimalDepth =
      j.at("expected_minimal_gates_at_minimal_depth").get<std::size_t>();
  test.expectedMinimalTwoQubitGates =
      j.at("expected_minimal_two_qubit_gates").get<std::size_t>();
  test.expectedMinimalGatesAtMinimalTwoQubitGates =
      j.at("expected_minimal_gates_at_minimal_two_qubit_gates")
          .get<std::size_t>();
}

namespace {
std::vector<TestConfiguration> getTests(const std::string& path) {
  std::ifstream input(path);
  nlohmann::json j;
  input >> j;
  return j;
}
} // namespace

CouplingMap parseEdges(const std::string& edgeString) {
  CouplingMap edges;
  std::stringstream ss(edgeString);
  std::string item;

  while (getline(ss, item, ';')) {
    item.erase(remove(item.begin(), item.end(), '{'), item.end());
    item.erase(remove(item.begin(), item.end(), '}'), item.end());
    const size_t pos = item.find(',');
    std::string first;
    std::string second;

    first = item.substr(0, pos);
    second = item.substr(pos + 1);

    const int u = stoi(first);
    const int v = stoi(second);

    // Insert the edge into the set
    edges.insert(Edge{u, v});
    std::cout << "Edge " << u << "," << v << " inserted to cm\n";
  }

  return edges;
}
class SynthesisTest : public ::testing::TestWithParam<TestConfiguration> {
protected:
  void SetUp() override {
    test = GetParam();

    if (!test.initialCircuit.empty()) {
      std::stringstream ss(test.initialCircuit);
      qc::QuantumComputation qc{};
      qc.import(ss, qc::Format::OpenQASM3);
      std::cout << "Initial circuit:\n" << qc;
      targetTableau = Tableau(qc);
      targetTableauWithDestabilizer =
          Tableau(qc, 0, std::numeric_limits<std::size_t>::max(), true);
      if (test.initialTableau.empty()) {
        initialTableau = Tableau(qc.getNqubits());
        initialTableauWithDestabilizer = Tableau(qc.getNqubits(), true);
        if (test.couplingMap.empty()) {
          synthesizer = CliffordSynthesizer(qc);
          synthesizerWithDestabilizer = CliffordSynthesizer(qc, true);
        } else {
          synthesizer = CliffordSynthesizer(qc, parseEdges(test.couplingMap));
          synthesizerWithDestabilizer =
              CliffordSynthesizer(qc, parseEdges(test.couplingMap), true);
        }
      } else {
        initialTableau = Tableau(test.initialTableau);
        std::cout << "Initial tableau:\n" << initialTableau;
        synthesizer = CliffordSynthesizer(initialTableau, qc);
      }
    } else {
      targetTableau = Tableau(test.targetTableau);
      if (test.initialTableau.empty()) {
        initialTableau = Tableau(targetTableau.getQubitCount());
        synthesizer = CliffordSynthesizer(targetTableau);
      } else {
        initialTableau = Tableau(test.initialTableau);
        std::cout << "Initial tableau:\n" << initialTableau;
        synthesizer = CliffordSynthesizer(initialTableau, targetTableau);
      }
    }
    std::cout << "Target tableau:\n" << targetTableau;

    config = Configuration();
    config.verbosity = plog::Severity::verbose;
    config.dumpIntermediateResults = true;
    config.useSymmetryBreaking = true;
  }

  void TearDown() override {
    std::cout << "Results:\n" << results << "\n";

    std::cout << "Mapping of Qubits:\n" << results.getMapping() << "\n";

    resultTableau = Tableau(synthesizer.getResultTableau());

    std::cout << "Result tableau:\n" << resultTableau;

    std::cout << "Target tableau:\n" << targetTableau;

    const std::vector<std::vector<bool>> p = results.getMappingVector();
    Tableau targetPrime = targetTableau.applyMapping(p);
    std::cout << "Target tableau with mapping:\n" << targetPrime;
    if (!targetPrime.hasDestabilizers()) {
      targetPrime.rref();
      resultTableau.rref();
      assert(targetPrime.equivalentUpToStabilizer(&resultTableau));
      std::cout << "Target tableau with mapping and Gauss:\n" << targetPrime;
      std::cout << "Result tableau with mapping and Gauss:\n" << resultTableau;
    } else {
      targetPrime =
          targetPrime.reverseMappingOnRows(p, targetPrime.getQubitCount());
      resultTableau =
          resultTableau.reverseMappingOnRows(p, targetPrime.getQubitCount());
      std::cout << "Result tableau with destab mapping reversed:\n"
                << resultTableau;
      std::cout << "Target tableau with destab mapping reversed:\n"
                << targetPrime;
    }
    EXPECT_EQ(resultTableau, targetPrime);

    const auto& resultCircuit = synthesizer.getResultCircuit();
    std::cout << "Resulting Circuit:\n" << resultCircuit;
    consistencyCheck(resultCircuit);
  }

  void consistencyCheck(const qc::QuantumComputation& qc) const {
    auto circuitTableau = initialTableau;
    for (const auto& gate : qc) {
      circuitTableau.applyGate(gate.get());
    }
    std::cout << "Circuit Tableau :\n" << circuitTableau;
    if (!circuitTableau.hasDestabilizers()) {
      circuitTableau.rref();
      std::cout << "Circuit Tableau with Gauss" << circuitTableau;
    }
    EXPECT_EQ(resultTableau, circuitTableau);
  }

  Tableau initialTableau;
  Tableau initialTableauWithDestabilizer;
  Tableau targetTableau;
  Tableau targetTableauWithDestabilizer;
  Configuration config;
  CliffordSynthesizer synthesizer;
  CliffordSynthesizer synthesizerWithDestabilizer;
  Results results;
  Results resultsWithDestabilizer;
  Tableau resultTableau;
  Tableau resultTableauWithDestabilizer;
  TestConfiguration test;
};

INSTANTIATE_TEST_SUITE_P(
    Tableaus, SynthesisTest, testing::ValuesIn(getTests("tableaus.json")),
    [](const testing::TestParamInfo<SynthesisTest::ParamType>& inf) {
      return inf.param.description;
    });

INSTANTIATE_TEST_SUITE_P(
    Circuits, SynthesisTest, testing::ValuesIn(getTests("circuits.json")),
    [](const testing::TestParamInfo<SynthesisTest::ParamType>& inf) {
      return inf.param.description;
    });

TEST_P(SynthesisTest, Gates) {
  config.target = TargetMetric::Gates;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), test.expectedMinimalGates);
}

TEST_P(SynthesisTest, GatesMaxSAT) {
  config.target = TargetMetric::Gates;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), test.expectedMinimalGates);
}

TEST_P(SynthesisTest, GatesLinearSearch) {
  config.target = TargetMetric::Gates;
  config.linearSearch = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), test.expectedMinimalGates);
}

TEST_P(SynthesisTest, Depth) {
  config.target = TargetMetric::Depth;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthMaxSAT) {
  config.target = TargetMetric::Depth;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthLinearSearch) {
  config.target = TargetMetric::Depth;
  config.linearSearch = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalGates) {
  config.target = TargetMetric::Depth;
  config.minimizeGatesAfterDepthOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
  EXPECT_EQ(results.getGates(), test.expectedMinimalGatesAtMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalTimeSteps) {
  config.target = TargetMetric::Depth;
  config.minimalTimesteps = test.expectedMinimalDepth;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalGatesMaxSAT) {
  config.target = TargetMetric::Depth;
  config.useMaxSAT = true;
  config.minimizeGatesAfterDepthOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
  EXPECT_EQ(results.getGates(), test.expectedMinimalGatesAtMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalGatesLinearSearch) {
  config.target = TargetMetric::Depth;
  config.linearSearch = true;
  config.minimizeGatesAfterDepthOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
  EXPECT_EQ(results.getGates(), test.expectedMinimalGatesAtMinimalDepth);
}

TEST_P(SynthesisTest, TwoQubitGates) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMaxSAT) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMinimalGates) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.minimizeGatesAfterTwoQubitGateOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
  EXPECT_EQ(results.getGates(),
            test.expectedMinimalGatesAtMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMinimalGatesMaxSAT) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.minimizeGatesAfterTwoQubitGateOptimization = true;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
  EXPECT_EQ(results.getGates(),
            test.expectedMinimalGatesAtMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TestDestabilizerGates) {
  if (!initialTableauWithDestabilizer.getTableau().empty()) {
    std::cout << "Testing with destabilizer\n";
    config.target = TargetMetric::Gates;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    synthesizerWithDestabilizer.synthesize(config);
    results = synthesizer.getResults();
    resultsWithDestabilizer = synthesizerWithDestabilizer.getResults();

    EXPECT_GE(resultsWithDestabilizer.getGates(), results.getGates());
  } else {
    std::cout << "Testing without destabilizer\n";
    config.target = TargetMetric::Gates;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    results = synthesizer.getResults();
  }
}

TEST_P(SynthesisTest, TestDestabilizerDepth) {
  if (!initialTableauWithDestabilizer.getTableau().empty()) {
    std::cout << "Testing with destabilizer\n";
    config.target = TargetMetric::Depth;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    synthesizerWithDestabilizer.synthesize(config);
    results = synthesizer.getResults();
    resultsWithDestabilizer = synthesizerWithDestabilizer.getResults();

    EXPECT_GE(resultsWithDestabilizer.getDepth(), results.getDepth());
  } else {
    std::cout << "Testing without destabilizer\n";
    config.target = TargetMetric::Gates;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    results = synthesizer.getResults();
  }
}

TEST_P(SynthesisTest, TestDestabilizerTwoQubitGates) {
  if (!initialTableauWithDestabilizer.getTableau().empty()) {
    std::cout << "Testing with destabilizer\n";
    config.target = TargetMetric::TwoQubitGates;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    synthesizerWithDestabilizer.synthesize(config);
    results = synthesizer.getResults();
    resultsWithDestabilizer = synthesizerWithDestabilizer.getResults();

    EXPECT_GE(resultsWithDestabilizer.getTwoQubitGates(),
              results.getTwoQubitGates());
  } else {
    std::cout << "Testing without destabilizer\n";
    config.target = TargetMetric::Gates;
    config.useMaxSAT = true;

    synthesizer.synthesize(config);
    results = synthesizer.getResults();
  }
}

TEST(HeuristicTest, basic) {
  auto config = Configuration();
  auto qc = qc::QuantumComputation(2);
  qc.h(0);
  qc.s(1);
  qc.h(0);
  qc.s(1);
  config.heuristic = true;
  config.splitSize = 1;
  config.target = TargetMetric::Depth;
  auto synth = CliffordSynthesizer(qc);
  synth.synthesize(config);
  EXPECT_EQ(synth.getResults().getDepth(), 2);
}

TEST(HeuristicTest, identity) {
  auto config = Configuration();
  auto qc = qc::QuantumComputation(2);
  qc.h(0);
  qc.s(1);
  qc.h(0);
  qc.sdg(1);
  config.heuristic = true;
  config.splitSize = 2;
  config.target = TargetMetric::Depth;
  auto synth = CliffordSynthesizer(qc);
  synth.synthesize(config);
  EXPECT_EQ(synth.getResults().getDepth(), 0);
}

TEST(HeuristicTest, threeLayers) {
  auto config = Configuration();
  auto qc = qc::QuantumComputation(2);
  qc.h(0);
  qc.h(1);
  qc.cx(0_pc, 1);
  qc.h(0);
  qc.h(1);
  config.heuristic = true;
  config.splitSize = 2;
  config.target = TargetMetric::Depth;
  auto synth = CliffordSynthesizer(qc);
  synth.synthesize(config);
  EXPECT_EQ(synth.getResults().getDepth(), 3);
}

TEST(HeuristicTest, fourLayers) {
  auto config = Configuration();
  auto qc = qc::QuantumComputation(1);
  qc.s(0);
  qc.s(0);
  qc.s(0);
  qc.s(0);
  config.heuristic = true;
  config.splitSize = 2;
  config.target = TargetMetric::Depth;
  auto synth = CliffordSynthesizer(qc);
  synth.synthesize(config);
  EXPECT_EQ(synth.getResults().getDepth(), 2);
}
} // namespace cs
