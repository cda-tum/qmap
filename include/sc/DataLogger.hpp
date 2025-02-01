//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Architecture.hpp"
#include "Definitions.hpp"
#include "MappingResults.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "utils.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

class DataLogger {
public:
  DataLogger(std::string path, Architecture& arch, qc::QuantumComputation qc)
      : dataLoggingPath(std::move(path)), architecture(&arch),
        nqubits(arch.getNqubits()), inputCircuit(std::move(qc)) {
    initLog();
    logArchitecture();
    logInputCircuit(inputCircuit);

    // collect registers
    auto combinedRegs = inputCircuit.getQuantumRegisters();
    for (const auto& reg : inputCircuit.getAncillaRegisters()) {
      combinedRegs.emplace(reg);
    }

    // build qubit index -> register map
    for (const auto& [_, reg] : combinedRegs) {
      const auto bound = reg.getStartIndex() + reg.getSize();
      for (qc::Qubit i = reg.getStartIndex(); i < bound; ++i) {
        qregs.try_emplace(i, reg, reg.toString(i));
      }
    }
    // build classical index -> register map
    for (const auto& [_, reg] : inputCircuit.getClassicalRegisters()) {
      const auto bound = reg.getStartIndex() + reg.getSize();
      for (qc::Bit i = reg.getStartIndex(); i < bound; ++i) {
        cregs.try_emplace(i, reg, reg.toString(i));
      }
    }
  }

  void initLog();
  void clearLog();
  void logArchitecture();
  void logSearchNode(std::size_t layer, std::size_t nodeId,
                     std::size_t parentId, double costFixed, double costHeur,
                     double lookaheadPenalty,
                     const std::vector<std::int16_t>& qubits, bool validMapping,
                     const std::vector<Exchange>& swaps, std::size_t depth);
  void logFinalizeLayer(
      std::size_t layer, const qc::CompoundOperation& ops,
      const std::vector<std::uint16_t>& singleQubitMultiplicity,
      const std::map<std::pair<std::uint16_t, std::uint16_t>,
                     std::pair<std::uint16_t, std::uint16_t>>&
          twoQubitMultiplicity,
      const std::vector<std::int16_t>& initialLayout, std::size_t finalNodeId,
      double finalCostFixed, double finalCostHeur, double finalLookaheadPenalty,
      const std::vector<std::int16_t>& finalLayout,
      const std::vector<Exchange>& finalSwaps, std::size_t finalSearchDepth);
  void splitLayer();
  void logMappingResult(MappingResults& result);
  void logInputCircuit(qc::QuantumComputation& qc) {
    if (deactivated) {
      return;
    }
    qc.dump(dataLoggingPath + "/input.qasm", qc::Format::OpenQASM3);
  };
  void logOutputCircuit(qc::QuantumComputation& qc) {
    if (deactivated) {
      return;
    }
    qc.dump(dataLoggingPath + "/output.qasm", qc::Format::OpenQASM3);
  }
  void close();

protected:
  std::string dataLoggingPath;
  Architecture* architecture;
  std::uint16_t nqubits;
  qc::QuantumComputation inputCircuit;
  qc::QubitIndexToRegisterMap qregs;
  qc::BitIndexToRegisterMap cregs;
  std::vector<std::ofstream> searchNodesLogFiles; // 1 per layer
  bool deactivated = false;

  void openNewLayer(std::size_t layer);
};
