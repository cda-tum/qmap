//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Definitions.hpp"
#include "cliffordsynthesis/CliffordSynthesizer.hpp"
#include "cliffordsynthesis/Configuration.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/TargetMetric.hpp"
#include "hybridmap/HybridNeutralAtomMapper.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomScheduler.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "na/NAComputation.hpp"
#include "na/NADefinitions.hpp"
#include "na/nasp/CodeGenerator.hpp"
#include "na/nasp/Solver.hpp"
#include "na/nasp/SolverFactory.hpp"
#include "sc/Architecture.hpp"
#include "sc/Mapper.hpp"
#include "sc/MappingResults.hpp"
#include "sc/configuration/AvailableArchitecture.hpp"
#include "sc/configuration/CommanderGrouping.hpp"
#include "sc/configuration/Configuration.hpp"
#include "sc/configuration/EarlyTermination.hpp"
#include "sc/configuration/Encoding.hpp"
#include "sc/configuration/Heuristic.hpp"
#include "sc/configuration/InitialLayout.hpp"
#include "sc/configuration/Layering.hpp"
#include "sc/configuration/LookaheadHeuristic.hpp"
#include "sc/configuration/Method.hpp"
#include "sc/configuration/SwapReduction.hpp"
#include "sc/exact/ExactMapper.hpp"
#include "sc/heuristic/HeuristicMapper.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <memory>
#include <nlohmann/json.hpp>
#include <plog/Severity.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11_json/pybind11_json.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

// c++ binding function
MappingResults map(const qc::QuantumComputation& circ, Architecture& arch,
                   Configuration& config) {
  if (config.useTeleportation) {
    config.teleportationQubits =
        std::min((arch.getNqubits() - circ.getNqubits()) & ~1U,
                 static_cast<std::size_t>(8));
  }

  std::unique_ptr<Mapper> mapper;
  try {
    if (config.method == Method::Heuristic) {
      mapper = std::make_unique<HeuristicMapper>(circ, arch);
    } else if (config.method == Method::Exact) {
      mapper = std::make_unique<ExactMapper>(circ, arch);
    }
  } catch (std::exception const& e) {
    std::stringstream ss{};
    ss << "Could not construct mapper: " << e.what();
    throw std::invalid_argument(ss.str());
  }

  try {
    mapper->map(config);
  } catch (std::exception const& e) {
    std::stringstream ss{};
    ss << "Error during mapping: " << e.what();
    throw std::invalid_argument(ss.str());
  }

  auto& results = mapper->getResults();

  std::stringstream qasm{};
  mapper->dumpResult(qasm, qc::Format::OpenQASM3);
  results.mappedCircuit = qasm.str();

  return results;
}

PYBIND11_MODULE(pyqmap, m, py::mod_gil_not_used()) {
  m.doc() = "pybind11 for the MQT QMAP quantum circuit mapping tool";

  // Pre-defined architecture available within QMAP
  py::enum_<AvailableArchitecture>(m, "Arch")
      .value("IBM_QX4", AvailableArchitecture::IbmQx4,
             "5 qubit, directed bow tie layout")
      .value("IBM_QX5", AvailableArchitecture::IbmQx5,
             "16 qubit, directed ladder layout")
      .value("IBMQ_Yorktown", AvailableArchitecture::IbmqYorktown,
             "5 qubit, undirected bow tie layout")
      .value("IBMQ_London", AvailableArchitecture::IbmqLondon,
             "5 qubit, undirected T-shape layout")
      .value("IBMQ_Bogota", AvailableArchitecture::IbmqBogota,
             "5 qubit, undirected linear chain layout")
      .value("IBMQ_Casablanca", AvailableArchitecture::IbmqCasablanca,
             "7 qubit, undirected H-shape layout")
      .value("IBMQ_Tokyo", AvailableArchitecture::IbmqTokyo,
             "20 qubit, undirected brick-like layout")
      .value("Rigetti_Agave", AvailableArchitecture::RigettiAgave,
             "8 qubit, undirected ring layout")
      .value("Rigetti_Aspen", AvailableArchitecture::RigettiAspen,
             "16 qubit, undirected dumbbell layout")
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> AvailableArchitecture {
        return architectureFromString(str);
      }));

  // Mapping methodology to use
  py::enum_<Method>(m, "Method")
      .value("heuristic", Method::Heuristic)
      .value("exact", Method::Exact)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Method {
        return methodFromString(str);
      }));

  // Initial layout strategy
  py::enum_<InitialLayout>(m, "InitialLayout")
      .value("identity", InitialLayout::Identity)
      .value("static", InitialLayout::Static)
      .value("dynamic", InitialLayout::Dynamic)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> InitialLayout {
        return initialLayoutFromString(str);
      }));

  // Heuristic function
  py::enum_<Heuristic>(m, "Heuristic")
      .value("gate_count_max_distance", Heuristic::GateCountMaxDistance)
      .value("gate_count_sum_distance", Heuristic::GateCountSumDistance)
      .value("gate_count_sum_distance_minus_shared_swaps",
             Heuristic::GateCountSumDistanceMinusSharedSwaps)
      .value("gate_count_max_distance_or_sum_distance_minus_shared_swaps",
             Heuristic::GateCountMaxDistanceOrSumDistanceMinusSharedSwaps)
      .value("fidelity_best_location", Heuristic::FidelityBestLocation)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Heuristic {
        return heuristicFromString(str);
      }));

  // Lookahead heuristic function
  py::enum_<LookaheadHeuristic>(m, "LookaheadHeuristic")
      .value("none", LookaheadHeuristic::None)
      .value("gate_count_max_distance",
             LookaheadHeuristic::GateCountMaxDistance)
      .value("gate_count_sum_distance",
             LookaheadHeuristic::GateCountSumDistance)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> LookaheadHeuristic {
        return lookaheadHeuristicFromString(str);
      }));

  // Gate clustering / layering strategy
  py::enum_<Layering>(m, "Layering")
      .value("individual_gates", Layering::IndividualGates)
      .value("disjoint_qubits", Layering::DisjointQubits)
      .value("odd_gates", Layering::OddGates)
      .value("qubit_triangle", Layering::QubitTriangle)
      .value("disjoint_2q_blocks", Layering::Disjoint2qBlocks)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Layering {
        return layeringFromString(str);
      }));

  // Early termination strategy in heuristic mapper
  py::enum_<EarlyTermination>(m, "EarlyTermination")
      .value("none", EarlyTermination::None)
      .value("expanded_nodes", EarlyTermination::ExpandedNodes)
      .value("expanded_nodes_after_first_solution",
             EarlyTermination::ExpandedNodesAfterFirstSolution)
      .value("expanded_nodes_after_current_optimal_solution",
             EarlyTermination::ExpandedNodesAfterCurrentOptimalSolution)
      .value("solution_nodes", EarlyTermination::SolutionNodes)
      .value("solution_nodes_after_current_optimal_solution",
             EarlyTermination::SolutionNodesAfterCurrentOptimalSolution)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> EarlyTermination {
        return earlyTerminationFromString(str);
      }));

  // Encoding settings for at-most-one and exactly-one constraints
  py::enum_<Encoding>(m, "Encoding")
      .value("naive", Encoding::Naive)
      .value("commander", Encoding::Commander)
      .value("bimander", Encoding::Bimander)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Encoding {
        return encodingFromString(str);
      }));

  // Grouping settings if using the commander encoding
  py::enum_<CommanderGrouping>(m, "CommanderGrouping")
      .value("fixed2", CommanderGrouping::Fixed2)
      .value("fixed3", CommanderGrouping::Fixed3)
      .value("halves", CommanderGrouping::Halves)
      .value("logarithm", CommanderGrouping::Logarithm)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> CommanderGrouping {
        return groupingFromString(str);
      }));

  // Strategy for reducing the number of permutations/swaps considered in front
  // of every gate
  py::enum_<SwapReduction>(m, "SwapReduction")
      .value("none", SwapReduction::None)
      .value("coupling_limit", SwapReduction::CouplingLimit)
      .value("custom", SwapReduction::Custom)
      .value("increasing", SwapReduction::Increasing)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> SwapReduction {
        return swapReductionFromString(str);
      }));

  // All configuration options for QMAP
  py::class_<Configuration>(
      m, "Configuration",
      "Configuration options for the MQT QMAP quantum circuit mapping tool")
      .def(py::init<>())
      .def_readwrite("method", &Configuration::method)
      .def_readwrite("heuristic", &Configuration::heuristic)
      .def_readwrite("verbose", &Configuration::verbose)
      .def_readwrite("debug", &Configuration::debug)
      .def_readwrite("data_logging_path", &Configuration::dataLoggingPath)
      .def_readwrite("layering", &Configuration::layering)
      .def_readwrite("automatic_layer_splits",
                     &Configuration::automaticLayerSplits)
      .def_readwrite("automatic_layer_splits_node_limit",
                     &Configuration::automaticLayerSplitsNodeLimit)
      .def_readwrite("early_termination", &Configuration::earlyTermination)
      .def_readwrite("early_termination_limit",
                     &Configuration::earlyTerminationLimit)
      .def_readwrite("initial_layout", &Configuration::initialLayout)
      .def_readwrite("iterative_bidirectional_routing",
                     &Configuration::iterativeBidirectionalRouting)
      .def_readwrite("iterative_bidirectional_routing_passes",
                     &Configuration::iterativeBidirectionalRoutingPasses)
      .def_readwrite("lookahead_heuristic", &Configuration::lookaheadHeuristic)
      .def_readwrite("lookaheads", &Configuration::nrLookaheads)
      .def_readwrite("first_lookahead_factor",
                     &Configuration::firstLookaheadFactor)
      .def_readwrite("lookahead_factor", &Configuration::lookaheadFactor)
      .def_readwrite("use_teleportation", &Configuration::useTeleportation)
      .def_readwrite("teleportation_qubits",
                     &Configuration::teleportationQubits)
      .def_readwrite("teleportation_seed", &Configuration::teleportationSeed)
      .def_readwrite("teleportation_fake", &Configuration::teleportationFake)
      .def_readwrite("timeout", &Configuration::timeout)
      .def_readwrite("encoding", &Configuration::encoding)
      .def_readwrite("commander_grouping", &Configuration::commanderGrouping)
      .def_readwrite("use_subsets", &Configuration::useSubsets)
      .def_readwrite("include_WCNF", &Configuration::includeWCNF)
      .def_readwrite("enable_limits", &Configuration::enableSwapLimits)
      .def_readwrite("swap_reduction", &Configuration::swapReduction)
      .def_readwrite("swap_limit", &Configuration::swapLimit)
      .def_readwrite("subgraph", &Configuration::subgraph)
      .def_readwrite("pre_mapping_optimizations",
                     &Configuration::preMappingOptimizations)
      .def_readwrite("post_mapping_optimizations",
                     &Configuration::postMappingOptimizations)
      .def_readwrite("add_measurements_to_mapped_circuit",
                     &Configuration::addMeasurementsToMappedCircuit)
      .def_readwrite("add_barriers_between_layers",
                     &Configuration::addBarriersBetweenLayers)
      .def("json", &Configuration::json)
      .def("__repr__", &Configuration::toString);

  // Results of the mapping process
  py::class_<MappingResults>(
      m, "MappingResults",
      "Results of the MQT QMAP quantum circuit mapping tool")
      .def(py::init<>())
      .def_readwrite("input", &MappingResults::input)
      .def_readwrite("output", &MappingResults::output)
      .def_readwrite("configuration", &MappingResults::config)
      .def_readwrite("time", &MappingResults::time)
      .def_readwrite("timeout", &MappingResults::timeout)
      .def_readwrite("mapped_circuit", &MappingResults::mappedCircuit)
      .def_readwrite("heuristic_benchmark", &MappingResults::heuristicBenchmark)
      .def_readwrite("layer_heuristic_benchmark",
                     &MappingResults::layerHeuristicBenchmark)
      .def_readwrite("wcnf", &MappingResults::wcnf)
      .def("json", &MappingResults::json)
      .def("csv", &MappingResults::csv)
      .def("__repr__", &MappingResults::toString);

  // Main class for storing circuit information
  py::class_<MappingResults::CircuitInfo>(m, "CircuitInfo",
                                          "Circuit information")
      .def(py::init<>())
      .def_readwrite("name", &MappingResults::CircuitInfo::name)
      .def_readwrite("qubits", &MappingResults::CircuitInfo::qubits)
      .def_readwrite("gates", &MappingResults::CircuitInfo::gates)
      .def_readwrite("single_qubit_gates",
                     &MappingResults::CircuitInfo::singleQubitGates)
      .def_readwrite("cnots", &MappingResults::CircuitInfo::cnots)
      .def_readwrite("layers", &MappingResults::CircuitInfo::layers)
      .def_readwrite("total_fidelity",
                     &MappingResults::CircuitInfo::totalFidelity)
      .def_readwrite("total_log_fidelity",
                     &MappingResults::CircuitInfo::totalLogFidelity)
      .def_readwrite("swaps", &MappingResults::CircuitInfo::swaps)
      .def_readwrite("direction_reverse",
                     &MappingResults::CircuitInfo::directionReverse)
      .def_readwrite("teleportations",
                     &MappingResults::CircuitInfo::teleportations);

  // Heuristic benchmark information
  py::class_<MappingResults::HeuristicBenchmarkInfo>(
      m, "HeuristicBenchmarkInfo", "Heuristic benchmark information")
      .def(py::init<>())
      .def_readwrite("expanded_nodes",
                     &MappingResults::HeuristicBenchmarkInfo::expandedNodes)
      .def_readwrite("generated_nodes",
                     &MappingResults::HeuristicBenchmarkInfo::generatedNodes)
      .def_readwrite("time_per_node",
                     &MappingResults::HeuristicBenchmarkInfo::secondsPerNode)
      .def_readwrite(
          "average_branching_factor",
          &MappingResults::HeuristicBenchmarkInfo::averageBranchingFactor)
      .def_readwrite(
          "effective_branching_factor",
          &MappingResults::HeuristicBenchmarkInfo::effectiveBranchingFactor)
      .def("json", &MappingResults::HeuristicBenchmarkInfo::json);

  // Heuristic benchmark information for individual layers
  py::class_<MappingResults::LayerHeuristicBenchmarkInfo>(
      m, "LayerHeuristicBenchmarkInfo", "Heuristic benchmark information")
      .def(py::init<>())
      .def_readwrite(
          "expanded_nodes",
          &MappingResults::LayerHeuristicBenchmarkInfo::expandedNodes)
      .def_readwrite(
          "generated_nodes",
          &MappingResults::LayerHeuristicBenchmarkInfo::generatedNodes)
      .def_readwrite("expanded_nodes_after_first_solution",
                     &MappingResults::LayerHeuristicBenchmarkInfo::
                         expandedNodesAfterFirstSolution)
      .def_readwrite("expanded_nodes_after_optimal_solution",
                     &MappingResults::LayerHeuristicBenchmarkInfo::
                         expandedNodesAfterOptimalSolution)
      .def_readwrite(
          "solution_nodes",
          &MappingResults::LayerHeuristicBenchmarkInfo::solutionNodes)
      .def_readwrite("solution_nodes_after_optimal_solution",
                     &MappingResults::LayerHeuristicBenchmarkInfo::
                         solutionNodesAfterOptimalSolution)
      .def_readwrite(
          "solution_depth",
          &MappingResults::LayerHeuristicBenchmarkInfo::solutionDepth)
      .def_readwrite(
          "time_per_node",
          &MappingResults::LayerHeuristicBenchmarkInfo::secondsPerNode)
      .def_readwrite(
          "average_branching_factor",
          &MappingResults::LayerHeuristicBenchmarkInfo::averageBranchingFactor)
      .def_readwrite("effective_branching_factor",
                     &MappingResults::LayerHeuristicBenchmarkInfo::
                         effectiveBranchingFactor)
      .def_readwrite(
          "early_termination",
          &MappingResults::LayerHeuristicBenchmarkInfo::earlyTermination)
      .def("json", &MappingResults::LayerHeuristicBenchmarkInfo::json);

  auto arch = py::class_<Architecture>(
      m, "Architecture", "Class representing device/backend information");
  auto properties = py::class_<Architecture::Properties>(
      arch, "Properties", "Class representing properties of an architecture");

  // Properties of an architecture (e.g. number of qubits, connectivity, error
  // rates, ...)
  properties.def(py::init<>())
      .def_property("name", &Architecture::Properties::getName,
                    &Architecture::Properties::setName)
      .def_property("num_qubits", &Architecture::Properties::getNqubits,
                    &Architecture::Properties::setNqubits)
      .def("get_single_qubit_error",
           &Architecture::Properties::getSingleQubitErrorRate, "qubit"_a,
           "operation"_a)
      .def("set_single_qubit_error",
           &Architecture::Properties::setSingleQubitErrorRate, "qubit"_a,
           "operation"_a, "error_rate"_a)
      .def("get_two_qubit_error",
           &Architecture::Properties::getTwoQubitErrorRate, "control"_a,
           "target"_a, "operation"_a = "cx")
      .def("set_two_qubit_error",
           &Architecture::Properties::setTwoQubitErrorRate, "control"_a,
           "target"_a, "error_rate"_a, "operation"_a = "cx")
      .def(
          "get_readout_error",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.readoutErrorRate.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_readout_error",
          [](Architecture::Properties& props, std::uint16_t qubit,
             double rate) { props.readoutErrorRate.set(qubit, rate); },
          "qubit"_a, "readout_error_rate"_a)
      .def(
          "get_t1",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.t1Time.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_t1",
          [](Architecture::Properties& props, std::uint16_t qubit, double t1) {
            props.t1Time.set(qubit, t1);
          },
          "qubit"_a, "t1"_a)
      .def(
          "get_t2",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.t2Time.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_t2",
          [](Architecture::Properties& props, std::uint16_t qubit, double t2) {
            props.t2Time.set(qubit, t2);
          },
          "qubit"_a, "t2"_a)
      .def(
          "get_frequency",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.qubitFrequency.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_frequency",
          [](Architecture::Properties& props, std::uint16_t qubit,
             double freq) { props.qubitFrequency.set(qubit, freq); },
          "qubit"_a, "qubit_frequency"_a)
      .def(
          "get_calibration_date",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.calibrationDate.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_calibration_date",
          [](Architecture::Properties& props, std::uint16_t qubit,
             const std::string& date) {
            props.calibrationDate.set(qubit, date);
          },
          "qubit"_a, "calibration_date"_a)
      .def("json", &Architecture::Properties::json,
           "Returns a JSON-style dictionary of all the information present in "
           "the :class:`.Properties`")
      .def("__repr__", &Architecture::Properties::toString,
           "Prints a JSON-formatted representation of all the information "
           "present in the :class:`.Properties`");

  // Interface to the QMAP internal architecture class
  arch.def(py::init<>())
      .def(py::init<std::uint16_t, const CouplingMap&>(), "num_qubits"_a,
           "coupling_map"_a)
      .def(py::init<std::uint16_t, const CouplingMap&,
                    const Architecture::Properties&>(),
           "num_qubits"_a, "coupling_map"_a, "properties"_a)
      .def_property("name", &Architecture::getName, &Architecture::setName)
      .def_property("num_qubits", &Architecture::getNqubits,
                    &Architecture::setNqubits)
      .def_property("coupling_map",
                    py::overload_cast<>(&Architecture::getCouplingMap),
                    &Architecture::setCouplingMap)
      .def_property("properties",
                    py::overload_cast<>(&Architecture::getProperties),
                    &Architecture::setProperties)
      .def("load_coupling_map",
           py::overload_cast<AvailableArchitecture>(
               &Architecture::loadCouplingMap),
           "available_architecture"_a)
      .def(
          "load_coupling_map",
          py::overload_cast<const std::string&>(&Architecture::loadCouplingMap),
          "coupling_map_file"_a)
      .def("load_properties",
           py::overload_cast<const Architecture::Properties&>(
               &Architecture::loadProperties),
           "properties"_a)
      .def("load_properties",
           py::overload_cast<const std::string&>(&Architecture::loadProperties),
           "properties"_a);

  // Main mapping function
  m.def("map", &map, "map a quantum circuit", "circ"_a, "arch"_a, "config"_a);

  // Target metric for the Clifford synthesizer
  py::enum_<cs::TargetMetric>(m, "TargetMetric")
      .value("gates", cs::TargetMetric::Gates, "Optimize gate count.")
      .value("two_qubit_gates", cs::TargetMetric::TwoQubitGates,
             "Optimize two-qubit gate count.")
      .value("depth", cs::TargetMetric::Depth, "Optimize circuit depth.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return cs::targetMetricFromString(name);
      }));
  py::implicitly_convertible<py::str, cs::TargetMetric>();

  py::enum_<plog::Severity>(m, "Verbosity")
      .value("none", plog::Severity::none, "No output.")
      .value("fatal", plog::Severity::fatal, "Only show fatal errors.")
      .value("error", plog::Severity::error, "Show errors.")
      .value("warning", plog::Severity::warning, "Show warnings.")
      .value("info", plog::Severity::info, "Show general information.")
      .value("debug", plog::Severity::debug,
             "Show additional debug information.")
      .value("verbose", plog::Severity::verbose, "Show all information.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return plog::severityFromString(name.c_str());
      }));
  py::implicitly_convertible<py::str, plog::Severity>();

  // Configuration for the synthesis
  py::class_<cs::Configuration>(
      m, "SynthesisConfiguration",
      "Configuration options for the MQT QMAP Clifford synthesis tool.")
      .def(py::init<>())
      .def_readwrite("initial_timestep_limit",
                     &cs::Configuration::initialTimestepLimit,
                     "Initial timestep limit for the Clifford synthesis. "
                     "Defaults to `0`, which implies that the initial timestep "
                     "limit is determined automatically.")
      .def_readwrite(
          "minimal_timesteps", &cs::Configuration::minimalTimesteps,
          "Minimal timestep considered for the Clifford synthesis. "
          "This option limits the lower bound of the interval in "
          "which the binary search method looks for solutions. "
          "Set this if you know a lower bound for the circuit depth. "
          "Defaults to `0`, which implies that no lower bound for depth "
          "is known.")
      .def_readwrite(
          "use_maxsat", &cs::Configuration::useMaxSAT,
          "Use MaxSAT to solve the synthesis problem or to really on the "
          "binary search scheme for finding the optimum. Defaults to `false`.")
      .def_readwrite("linear_search", &cs::Configuration::linearSearch,
                     "Use liner search instead of binary search "
                     "scheme for finding the optimum. Defaults to `false`.")
      .def_readwrite(
          "target_metric", &cs::Configuration::target,
          "Target metric for the Clifford synthesis. Defaults to `gates`.")
      .def_readwrite("use_symmetry_breaking",
                     &cs::Configuration::useSymmetryBreaking,
                     "Use symmetry breaking clauses to speed up the synthesis "
                     "process. Defaults to `true`.")
      .def_readwrite("dump_intermediate_results",
                     &cs::Configuration::dumpIntermediateResults,
                     "Dump intermediate results of the synthesis process. "
                     "Defaults to `false`.")
      .def_readwrite("intermediate_results_path",
                     &cs::Configuration::intermediateResultsPath,
                     "Path to the directory where intermediate results should "
                     "be dumped. Defaults to `./`. The path needs to include a "
                     "path separator at the end.")
      .def_readwrite(
          "verbosity", &cs::Configuration::verbosity,
          "Verbosity level for the synthesis process. Defaults to 'warning'.")
      .def_readwrite("solver_parameters", &cs::Configuration::solverParameters,
                     "Parameters to be passed to Z3 as dict[str, bool | int | "
                     "float | str]")
      .def_readwrite(
          "minimize_gates_after_depth_optimization",
          &cs::Configuration::minimizeGatesAfterDepthOptimization,
          "Depth optimization might produce a circuit with more gates than "
          "necessary. This option enables an additional run of the synthesizer "
          "to minimize the overall number of gates. Defaults to `false`.")
      .def_readwrite(
          "try_higher_gate_limit_for_two_qubit_gate_optimization",
          &cs::Configuration::tryHigherGateLimitForTwoQubitGateOptimization,
          "When optimizing two-qubit gates, the synthesizer might fail "
          "to find an optimal solution for a certain timestep limit, but there "
          "might be a better solution for some higher timestep limit. This "
          "option enables an additional run of the synthesizer with a higher "
          "gate limit. Defaults to `false`.")
      .def_readwrite("gate_limit_factor", &cs::Configuration::gateLimitFactor,
                     "Factor by which the gate limit is increased when "
                     "trying to find a better solution for the two-qubit "
                     "gate optimization. Defaults to `1.1`.")
      .def_readwrite(
          "minimize_gates_after_two_qubit_gate_optimization",
          &cs::Configuration::minimizeGatesAfterTwoQubitGateOptimization,
          "Two-qubit gate optimization might produce a circuit "
          "with more gates than necessary. This option enables "
          "an additional run of the synthesizer to minimize the "
          "overall number of gates. Defaults to `false`.")
      .def_readwrite("heuristic", &cs::Configuration::heuristic,
                     "Use heuristic to synthesize the circuit. "
                     "This method synthesizes shallow intermediate circuits "
                     "and combines them. Defaults to `false`.")
      .def_readwrite("split_size", &cs::Configuration::splitSize,
                     "Size of subcircuits used in heuristic. "
                     "Defaults to `5`.")
      .def_readwrite(
          "n_threads_heuristic", &cs::Configuration::nThreadsHeuristic,
          "Maximum number of threads used for the heuristic optimizer. "
          "Defaults to the number of available threads on the system.")
      .def("json", &cs::Configuration::json,
           "Returns a JSON-style dictionary of all the information present in "
           "the :class:`.Configuration`")
      .def(
          "__repr__",
          [](const cs::Configuration& config) { return config.json().dump(2); },
          "Prints a JSON-formatted representation of all the information "
          "present in the :class:`.Configuration`");

  // Results of the synthesis
  py::class_<cs::Results>(m, "SynthesisResults",
                          "Results of the MQT QMAP Clifford synthesis tool.")
      .def(py::init<>())
      .def_property_readonly("gates", &cs::Results::getGates,
                             "Returns the number of gates in the circuit.")
      .def_property_readonly("single_qubit_gates",
                             &cs::Results::getSingleQubitGates,
                             "Returns the number of single-qubit gates in the "
                             "synthesized circuit.")
      .def_property_readonly("two_qubit_gates", &cs::Results::getTwoQubitGates,
                             "Returns the number of two-qubit gates in the "
                             "synthesized circuit.")
      .def_property_readonly("depth", &cs::Results::getDepth,
                             "Returns the depth of the synthesized circuit.")
      .def_property_readonly("runtime", &cs::Results::getRuntime,
                             "Returns the runtime of the synthesis in seconds.")
      .def_property_readonly("solver_calls", &cs::Results::getSolverCalls,
                             "Returns the number of calls to the SAT solver.")
      .def_property_readonly(
          "circuit", &cs::Results::getResultCircuit,
          "Returns the synthesized circuit as a qasm string.")
      .def_property_readonly("tableau", &cs::Results::getResultTableau,
                             "Returns a string representation of the "
                             "synthesized circuit's tableau.")
      .def("sat", &cs::Results::sat,
           "Returns `true` if the synthesis was successful.")
      .def("unsat", &cs::Results::unsat,
           "Returns `true` if the synthesis was unsuccessful.");

  auto tableau = py::class_<cs::Tableau>(
      m, "Tableau", "A class for representing stabilizer tableaus.");
  tableau.def(py::init<std::size_t, bool>(), "n"_a,
              "include_destabilizers"_a = false,
              "Creates a tableau for an n-qubit Clifford.");
  tableau.def(
      py::init<const std::string&>(), "tableau"_a,
      "Constructs a tableau from a string description. This can either be a "
      "semicolon separated binary matrix or a list of Pauli strings.");
  tableau.def(
      py::init<const std::string&, const std::string&>(), "stabilizers"_a,
      "destabilizers"_a,
      "Constructs a tableau from two lists of Pauli strings, the Stabilizers"
      "and Destabilizers.");

  auto synthesizer = py::class_<cs::CliffordSynthesizer>(
      m, "CliffordSynthesizer", "A class for synthesizing Clifford circuits.");

  synthesizer.def(py::init<cs::Tableau, cs::Tableau>(), "initial_tableau"_a,
                  "target_tableau"_a,
                  "Constructs a synthesizer for two tableaus representing the "
                  "initial and target state.");
  synthesizer.def(py::init<cs::Tableau>(), "target_tableau"_a,
                  "Constructs a synthesizer for a tableau representing the "
                  "target state.");
  synthesizer.def(py::init<qc::QuantumComputation&, bool>(), "qc"_a,
                  "use_destabilizers"_a
                  "Constructs a synthesizer for a quantum computation "
                  "representing the target state.");
  synthesizer.def(
      py::init<cs::Tableau, qc::QuantumComputation&>(), "initial_tableau"_a,
      "qc"_a,
      "Constructs a synthesizer for a quantum computation representing the "
      "target state that starts in an initial state represented by a tableau.");
  synthesizer.def("synthesize", &cs::CliffordSynthesizer::synthesize,
                  "config"_a = cs::Configuration(),
                  "Runs the synthesis with the given configuration.");
  synthesizer.def_property_readonly("results",
                                    &cs::CliffordSynthesizer::getResults,
                                    "Returns the results of the synthesis.");

  // Neutral Atom Hybrid Mapper
  py::enum_<na::InitialCoordinateMapping>(
      m, "InitialCoordinateMapping",
      "Initial mapping between hardware qubits hardware coordinates.")
      .value("trivial", na::InitialCoordinateMapping::Trivial,
             "Trivial identity mapping.")
      .value("random", na::InitialCoordinateMapping::Random, "Random mapping.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return na::initialCoordinateMappingFromString(name);
      }));
  py::enum_<na::InitialMapping>(
      m, "InitialCircuitMapping",
      "Initial mapping between circuit qubits and hardware qubits.")
      .value("identity", na::InitialMapping::Identity, "Identity mapping.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return na::initialMappingFromString(name);
      }));

  py::class_<na::MapperParameters>(
      m, "HybridMapperParameters",
      "Parameters for the Neutral Atom Hybrid Mapper.")
      .def(py::init<>())
      .def(py::init<double, double, double, double, double, double>(),
           "lookahead_weight_swaps"_a = 0.1, "lookahead_weight_moves"_a = 0.1,
           "decay"_a = 0.1, "shuttling_time_weight"_a = 1, "gate_weight"_a = 1,
           "shuttling_weight"_a = 1)
      .def_readwrite("lookahead_weight_swaps",
                     &na::MapperParameters::lookaheadWeightSwaps,
                     "Weight for the lookahead for the SWAP gates. 0 means no "
                     "lookahead is considered.")
      .def_readwrite("lookahead_weight_moves",
                     &na::MapperParameters::lookaheadWeightMoves,
                     "Weight for the lookahead for the MOVE gates. 0 means no "
                     "lookahead is considered.")
      .def_readwrite("decay", &na::MapperParameters::decay,
                     "Decay factor for the blocking constraint to avoid SWAPs "
                     "that block each other. 0 means no decay.")
      .def_readwrite(
          "shuttling_time_weight", &na::MapperParameters::shuttlingTimeWeight,
          "Weight how much the shuttling Times should be considered.")
      .def_readwrite(
          "gate_weight", &na::MapperParameters::gateWeight,
          "Weight for the SWAP gates. Higher means mapper will prefer SWAP "
          "gates over shuttling. 0 means only shuttling is used.")
      .def_readwrite(
          "shuttling_weight", &na::MapperParameters::shuttlingWeight,
          "Weight for the shuttling. Higher means mapper will prefer "
          "shuttling over SWAP gates. 0 means only SWAP gates are "
          "used.")
      .def_readwrite(
          "seed", &na::MapperParameters::seed,
          "Seed for the random number generator. 0 means random seed.")
      .def_readwrite("verbose", &na::MapperParameters::verbose,
                     "Print additional information during the mapping process.")
      .def_readwrite("initial_mapping", &na::MapperParameters::initialMapping,
                     "Initial mapping between circuit qubits and hardware "
                     "qubits.");

  py::class_<na::NeutralAtomArchitecture>(m, "NeutralAtomHybridArchitecture")
      .def(py::init<const std::string&>(), "filename"_a)
      .def("load_json", &na::NeutralAtomArchitecture::loadJson,
           "json_filename"_a)
      .def_readwrite("name", &na::NeutralAtomArchitecture::name,
                     "Name of the "
                     "architecture")
      .def_property_readonly(
          "nrows", &na::NeutralAtomArchitecture::getNrows,
          "Number of rows in a rectangular grid SLM arrangement")
      .def_property_readonly(
          "ncolumns", &na::NeutralAtomArchitecture::getNcolumns,
          "Number of columns in a rectangular grid SLM arrangement")
      .def_property_readonly(
          "npositions", &na::NeutralAtomArchitecture::getNpositions,
          "Total number of positions in a rectangular grid SLM arrangement")
      .def_property_readonly(
          "naods", &na::NeutralAtomArchitecture::getNAods,
          "Number of independent 2D acousto-optic deflectors")
      .def_property_readonly("naod_coordinates",
                             &na::NeutralAtomArchitecture::getNAodCoordinates,
                             "Maximal number of AOD rows/columns (NOT USED)")
      .def_property_readonly("nqubits",
                             &na::NeutralAtomArchitecture::getNqubits,
                             "Number of atoms in the neutral atom quantum "
                             "computer that can be used as qubits")
      .def_property_readonly(
          "inter_qubit_distance",
          &na::NeutralAtomArchitecture::getInterQubitDistance,
          "Distance "
          "between "
          "SLM traps in "
          "micrometers")
      .def_property_readonly("interaction_radius",
                             &na::NeutralAtomArchitecture::getInteractionRadius,
                             "Interaction radius in inter-qubit distances")
      .def_property_readonly("blocking_factor",
                             &na::NeutralAtomArchitecture::getBlockingFactor,
                             "Blocking factor for parallel Rydberg gates")
      .def_property_readonly(
          "naod_intermediate_levels",
          &na::NeutralAtomArchitecture::getNAodIntermediateLevels,
          "Number of possible AOD positions between two SLM traps")
      .def_property_readonly("decoherence_time",
                             &na::NeutralAtomArchitecture::getDecoherenceTime,
                             "Decoherence time in microseconds")
      .def("compute_swap_distance",
           static_cast<int (na::NeutralAtomArchitecture::*)(
               std::uint32_t, std::uint32_t) const>(
               &na::NeutralAtomArchitecture::getSwapDistance),
           "Number of SWAP gates required between two positions",
           py::arg("idx1"), py::arg("idx2"))
      .def("get_gate_time", &na::NeutralAtomArchitecture::getGateTime,
           "Execution time of certain gate in microseconds", "s"_a)
      .def("get_gate_average_fidelity",
           &na::NeutralAtomArchitecture::getGateAverageFidelity,
           "Average gate fidelity from [0,1]", "s"_a)
      .def("get_nearby_coordinates",
           &na::NeutralAtomArchitecture::getNearbyCoordinates,
           "Positions that are within the interaction radius of the passed "
           "position",
           "idx"_a)
      .def("get_animation_csv", &na::NeutralAtomArchitecture::getAnimationCsv,
           "Returns string representation of the architecture used for "
           "animation")
      .def("save_animation_csv", &na::NeutralAtomArchitecture::saveAnimationCsv,
           "filename"_a, "Saves the animation csv string to a file");

  py::class_<na::NeutralAtomMapper>(
      m, "HybridNAMapper",
      "Neutral Atom Hybrid Mapper that can use both SWAP gates and AOD "
      "movements to map a quantum circuit to a neutral atom quantum computer.")
      .def(py::init<const na::NeutralAtomArchitecture&, na::MapperParameters>(),
           "Create Hybrid NA Mapper with mapper parameters",
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), "arch"_a,
           "params"_a = na::MapperParameters())
      .def("set_parameters", &na::NeutralAtomMapper::setParameters,
           "Set the parameters for the Hybrid NA Mapper", "params"_a)
      .def(
          "get_init_hw_pos", &na::NeutralAtomMapper::getInitHwPos,
          "Get the initial hardware positions, required to create an animation")
      .def("map", &na::NeutralAtomMapper::mapAndConvert,
           "Map a quantum circuit to the neutral atom quantum computer",
           "circ"_a, "initial_mapping"_a = na::InitialMapping::Identity,
           "verbose"_a = false)
      .def(
          "map_qasm_file",
          [](na::NeutralAtomMapper& mapper, const std::string& filename,
             const na::InitialMapping initialMapping) {
            qc::QuantumComputation qc(filename);
            mapper.map(qc, initialMapping);
          },
          "Map a quantum circuit to the neutral atom quantum computer",
          "filename"_a, "initial_mapping"_a = na::InitialMapping::Identity)
      .def("get_mapped_qc", &na::NeutralAtomMapper::getMappedQc,
           "Returns the mapped circuit as an extended qasm2 string")
      .def("save_mapped_qc", &na::NeutralAtomMapper::saveMappedQc,
           "Saves the mapped circuit as an extended qasm2 string to a file",
           "filename"_a)
      .def("get_mapped_qc_aod", &na::NeutralAtomMapper::getMappedQcAOD,
           "Returns the mapped circuit as an extended qasm2 string with native "
           "AOD movements")
      .def("save_mapped_qc_aod", &na::NeutralAtomMapper::saveMappedQcAOD,
           "Saves the mapped circuit as an extended qasm2 string with native "
           "AOD movements to a file",
           "filename"_a)
      .def(
          "schedule",
          [](na::NeutralAtomMapper& mapper, const bool verbose,
             const bool createAnimationCsv, const double shuttlingSpeedFactor) {
            auto results = mapper.schedule(verbose, createAnimationCsv,
                                           shuttlingSpeedFactor);
            return results.toMap();
          },
          "Schedule the mapped circuit", "verbose"_a = false,
          "create_animation_csv"_a = false, "shuttling_speed_factor"_a = 1.0)
      .def("get_animation_csv", &na::NeutralAtomMapper::getAnimationCsv,
           "Returns the animation csv string")
      .def("save_animation_csv", &na::NeutralAtomMapper::saveAnimationCsv,
           "Saves the animation csv string to a file", "filename"_a);

  // Neutral Atom State Preparation
  py::class_<na::NASolver>(m, "NAStatePreparationSolver", R"(
The neutral atom state preparation solver generates an optimal sequence of
neutral atom operations for a given state preparation circuit.
)")
      .def(py::init<uint16_t, uint16_t, uint16_t, uint16_t, uint16_t, uint16_t,
                    uint16_t, uint16_t, uint16_t, uint16_t>(),
           "max_x"_a, "max_y"_a, "max_c"_a, "max_r"_a, "max_h_offset"_a,
           "max_v_offset"_a, "max_h_dist"_a, "max_v_dist"_a,
           "min_entangling_y"_a, "max_entangling_y"_a, R"(
Create a solver instance for the neutral atom state preparation problem.

The solver is based on a 2D grid abstraction of the neutral atom quantum
computer. The 2D plane is divided into interaction sites. Each interaction site
is denoted by abstract x- and y-coordinates. The parameter `max_x` specifies the
maximum x-coordinate, and `max_y` specifies the maximum y-coordinate. In the
center of an interaction site, sits an SLM trap. Around that trap there are
several possible discrete AOD positions arranged as a grid. The specific
position of an atom within an interaction site is determined by the x- and
y-offset from the SLM trap, i.e., those can be positive and negative. The
maximum absolute value of the x- and y-offset is specified by `max_h_offset` and
`max_v_offset`, respectively. Then, the parameter `max_c` specifies the maximum
number of AOD columns, and `max_r` specifies the maximum number of AOD rows.
Finally, in order to interact during a Rydberg stage, atoms must be located
within a certain distance. The maximum horizontal and vertical distance between
two atoms is specified by `max_h_dist` and `max_v_dist`, respectively. The
parameter `min_entangling_y` specifies the minimum y-coordinate for entangling
operations, and `max_entangling_y` specifies the maximum y-coordinate for
entangling operations. Hence, y-coordinates outside of this range are located in
the storage zone.

.. note::
    The solver can only handle a single storage zone below the entangling zone,
    i.e., in this case `min_entangling_y` must be zero and `max_entangling_y`
    must be less than `max_y`.

:param max_x: is the maximum discrete x-coordinate of the interaction sites
:param max_y: is the maximum discrete y-coordinate of the interaction sites
:param max_c: is the maximum number of AOD columns
:param max_r: is the maximum number of AOD rows
:param max_h_offset: is the maximum horizontal offset of the atoms
:param max_v_offset: is the maximum vertical offset of the atoms
:param max_h_dist: is the maximum horizontal distance between two atoms
:param max_v_dist: is the maximum vertical distance between two atoms
:param min_entangling_y: is the minimum y-coordinate for entangling operations
:param max_entangling_y: is the maximum y-coordinate for entangling operations
:raises ValueError: if one of the parameters is invalid, e.g., is a negative
value
)")
      .def("solve", &na::NASolver::solve, "ops"_a, "num_qubits"_a,
           "num_stages"_a, "num_transfers"_a, "mind_ops_order"_a,
           "shield_idle_qubits"_a, R"(
Solve the neutral atom state preparation problem.

The solver generates an optimal sequence of neutral atom operations for a given
state preparation circuit. The circuit is given as a list of operations, where
each operation is a pair of qubits. The sequence is divided into stages. Each
stage is either a Rydberg stage or a transfer stage. In a Rydberg stage,
adjacent qubits in the entangling zone undergo an entangling gate. In a transfer
stage, atoms can be stored from AOD into SLM traps and loaded from SLM traps
into AOD. At the end of each stage, the atoms are shuttled to their next
position. The number of stages is specified by `num_stages`. The number of
transfers is fixed by `num_transfers` if give. If this parameter is not
specified, then the solver will determine the optimal number of transfers. The
parameter `mind_ops_order` specifies whether the order of the operations in the
circuit should be preserved. The parameter `shield_idle_qubits` specifies
whether idle qubits should be shielded from the entangling operations.

.. note::
    To retrieve the list of qubit pairs from a quantum circuit, use the function
    :func:`get_ops_for_solver`.

.. note::
    The returned solver's result can either directly exported to the YAML format
    by calling the method :func:`yaml` on the result object or the result object
    can be passed to the function :func:`generate_code` to generate code
    consisting of neutral atom operations.

:param ops: is the list of operations in the circuit
:param num_qubits: is the number of qubits in the circuit
:param num_stages: is the number of stages in the sequence
:param num_transfers: (optional) is the number of transfers in the sequence
:param mind_ops_order: is True if the order of the operations should be
  preserved
:param shield_idle_qubits: is True if idle qubits should be shielded
:returns: the result of the solver
:raises ValueError: if one of the numeral parameters is invalid, e.g., is a
  negative value
)");

  py::class_<na::NASolver::Result>(
      m, "NAStatePreparationSolver.Result",
      "Neutral Atom State Preparation Solver Result")
      .def(py::init<>(), "Create a result object")
      .def(
          "yaml",
          [](const na::NASolver::Result& result, const bool compact) {
            return result.yaml(0, compact);
          },
          "compact"_a = true, R"(
Returns the result as a YAML string.

:param compact: if True, the YAML string contains JSON elements to make it more
compact
:returns: the result as a YAML string
)");

  m.def(
      "generate_code",
      [](const qc::QuantumComputation& qc, const na::NASolver::Result& result,
         const uint16_t minAtomDist, const uint16_t noInteractionRadius,
         const uint16_t zoneDist) {
        return na::CodeGenerator::generate(qc, result, minAtomDist,
                                           noInteractionRadius, zoneDist)
            .toString();
      },
      "qc"_a, "result"_a, "min_atom_dist"_a = 1, "no_interaction_radius"_a = 10,
      "zone_dist"_a = 24, R"(
Generate code for the given circuit using the solver's result. Some parameters
of the abstraction from the 2D grid used for the solver must be provided again.

:param qc: is the quantum circuit
:param result: is the result of the solver
:param min_atom_dist: is the minimum distance between atoms
:param no_interaction_radius: is the radius around an atom where no other atom
  can be placed during an entangling operation that should not interact with the
  atom
:param zone_dist: is the distance between zones, i.e., the minimal distance
  between two atoms in different zones
:raises ValueError: if one of the numeral parameters is invalid, e.g., is a
  negative value)");

  m.def(
      "get_ops_for_solver",
      [](const qc::QuantumComputation& qc, const std::string& operationType,
         const uint64_t numControls, const bool quiet) {
        auto opTypeLowerStr = operationType;
        std::transform(opTypeLowerStr.begin(), opTypeLowerStr.end(),
                       opTypeLowerStr.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        const auto fullOpType =
            na::FullOpType{qc::opTypeFromString(operationType), numControls};
        return na::SolverFactory::getOpsForSolver(qc, fullOpType, quiet);
      },
      "qc"_a, "operation_type"_a = "Z", "num_operands"_a = 1, "quiet"_a = true,
      R"(
Extract entangling operations as list of qubit pairs from the circuit.

.. warning::
    This function can only extract qubit pairs of two-qubit operations.
    I.e., the operands of the operation plus the controls must be equal to two.

:param qc: is the quantum circuit
:param operation_type: is the type of operation to extract, e.g., "z" for CZ
  gates
:param num_controls: is the number of controls the operation acts on, e.g., 1
  for CZ gates
:param quiet: if True, suppresses warning when the circuit contains operations
  other than the specified operation type
:return: list of qubit pairs
:raises ValueError: if the circuit contains operations other than the specified
  operation type and quiet is False
:raises ValueError: if the operation has more than two operands including
  controls
)");
}
