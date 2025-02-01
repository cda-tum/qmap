//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "sc/Mapper.hpp"

#include "Definitions.hpp"
#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "sc/Architecture.hpp"
#include "sc/configuration/Layering.hpp"
#include "sc/utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <utility>
#include <vector>

void Mapper::initResults() {
  countGates(qc, results.input);
  results.input.name = qc.getName();
  results.input.qubits = static_cast<std::uint16_t>(qc.getNqubits());
  results.architecture = architecture->getName();
  results.output.name = qc.getName() + "_mapped";
  results.output.qubits = architecture->getNqubits();
  results.output.gates = std::numeric_limits<std::size_t>::max();
  qcMapped.addQubitRegister(architecture->getNqubits());
}

Mapper::Mapper(qc::QuantumComputation quantumComputation, Architecture& arch)
    : qc(std::move(quantumComputation)), architecture(&arch) {
  qubits.resize(architecture->getNqubits(), DEFAULT_POSITION);
  locations.resize(architecture->getNqubits(), DEFAULT_POSITION);

  // strip away qubits that are not used in the circuit
  qc.stripIdleQubits(true);
  // strip away final measurement gates
  qc::CircuitOptimizer::removeFinalMeasurements(qc);
}

void Mapper::processDisjointQubitLayer(
    std::vector<std::optional<std::size_t>>& lastLayer,
    const std::optional<std::uint16_t>& control, const std::uint16_t target,
    qc::Operation* gate) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target) + 1;
    }
    lastLayer.at(target) = layer;
  } else {
    if (!lastLayer.at(*control).has_value() &&
        !lastLayer.at(target).has_value()) {
      layer = 0;
    } else if (!lastLayer.at(*control).has_value()) {
      layer = *lastLayer.at(target) + 1;
    } else if (!lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(*control) + 1;
    } else {
      layer = std::max(*lastLayer.at(*control), *lastLayer.at(target)) + 1;
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target) = layer;
  }

  if (layers.size() <= layer) {
    layers.emplace_back();
  }
  if (control.has_value()) {
    layers.at(layer).emplace_back(*control, target, gate);
  } else {
    layers.at(layer).emplace_back(-1, target, gate);
  }
}

void Mapper::processDisjoint2qBlockLayer(
    std::vector<std::optional<std::size_t>>& lastLayer,
    const std::optional<std::uint16_t>& control, const std::uint16_t target,
    qc::Operation* gate) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    // single qubit gates can always be added to the last 2Q block and should
    // not affect placings of future 2Q blocks
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target);
    }
  } else {
    if (!lastLayer.at(*control).has_value() &&
        !lastLayer.at(target).has_value()) {
      layer = 0;
    } else if (!lastLayer.at(*control).has_value()) {
      layer = *lastLayer.at(target) + 1;
    } else if (!lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(*control) + 1;
    } else {
      layer = std::max(*lastLayer.at(*control), *lastLayer.at(target)) + 1;

      if (*lastLayer.at(*control) == *lastLayer.at(target)) {
        for (auto& g : layers.at(layer - 1)) {
          if ((g.control == *control && g.target == target) ||
              (g.control == target && g.target == *control)) {
            // if last layer contained gate with equivalent qubit set, use that
            // layer
            --layer;
            break;
          }
        }
      }
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target) = layer;
  }

  if (layers.size() <= layer) {
    layers.emplace_back();
  }
  if (control.has_value()) {
    layers.at(layer).emplace_back(*control, target, gate);
  } else {
    layers.at(layer).emplace_back(-1, target, gate);
  }
}

void Mapper::createLayers() {
  const auto& config = results.config;
  auto lastLayer = std::vector<std::optional<std::size_t>>(
      architecture->getNqubits(), std::nullopt);

  auto qubitsInLayer = std::set<std::uint16_t>{};

  bool even = true;
  for (auto& gate : qc) {
    // skip over barrier instructions
    if (gate->getType() == qc::Barrier || gate->getType() == qc::Measure) {
      continue;
    }

    if (!gate->isUnitary()) {
      throw QMAPException(
          "Mapping not possible: circuit contains non-unitary operation: " +
          std::string(gate->getName()));
    }

    if (gate->getControls().size() > 1 || gate->getTargets().size() > 1) {
      throw QMAPException("Circuit contains gates with more than one control. "
                          "Please make sure that the input circuit's gates are "
                          "decomposed to the appropriate gate set!");
    }

    const bool singleQubit = !gate->isControlled();
    std::optional<std::uint16_t> control = std::nullopt;
    if (!singleQubit) {
      control = static_cast<std::uint16_t>(
          qc.initialLayout.at((*gate->getControls().begin()).qubit));
    }
    const auto target = static_cast<std::uint16_t>(
        qc.initialLayout.at(gate->getTargets().at(0)));

    // methods of layering described in
    // https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf
    switch (config.layering) {
    case Layering::IndividualGates:
      // each gate is put in a new layer
      layers.emplace_back();
      if (control.has_value()) {
        layers.back().emplace_back(*control, target, gate.get());
      } else {
        layers.back().emplace_back(-1, target, gate.get());
      }
      break;
    case Layering::DisjointQubits:
      processDisjointQubitLayer(lastLayer, control, target, gate.get());
      break;
    case Layering::Disjoint2qBlocks:
      processDisjoint2qBlockLayer(lastLayer, control, target, gate.get());
      break;
    case Layering::OddGates:
      // every other gate is put in a new layer
      if (even) {
        layers.emplace_back();
        if (control.has_value()) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.back().emplace_back(-1, target, gate.get());
        }
      } else {
        if (control.has_value()) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.back().emplace_back(-1, target, gate.get());
        }
      }
      even = !even;
      break;
    case Layering::QubitTriangle:
      if (layers.empty()) {
        layers.emplace_back();
      }

      if (singleQubit) {
        // single qubit gates can be added in any layer
        layers.back().emplace_back(-1, target, gate.get());
      } else {
        qubitsInLayer.insert(*control);
        qubitsInLayer.insert(target);

        if (qubitsInLayer.size() <= 3) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.emplace_back();
          layers.back().emplace_back(*control, target, gate.get());
          qubitsInLayer.clear();
          qubitsInLayer.insert(*control);
          qubitsInLayer.insert(target);
        }
      }
      break;
    }
  }
  results.input.layers = layers.size();

  // compute qubit gate multiplicities
  singleQubitMultiplicities = std::vector<SingleQubitMultiplicity>(
      layers.size(), SingleQubitMultiplicity(architecture->getNqubits(), 0));
  twoQubitMultiplicities =
      std::vector<TwoQubitMultiplicity>(layers.size(), TwoQubitMultiplicity{});
  activeQubits = std::vector<std::set<std::uint16_t>>(
      layers.size(), std::set<std::uint16_t>{});
  activeQubits1QGates = std::vector<std::set<std::uint16_t>>(
      layers.size(), std::set<std::uint16_t>{});
  activeQubits2QGates = std::vector<std::set<std::uint16_t>>(
      layers.size(), std::set<std::uint16_t>{});

  for (std::size_t i = 0; i < layers.size(); ++i) {
    for (const auto& gate : layers[i]) {
      if (gate.singleQubit()) {
        activeQubits[i].emplace(gate.target);
        activeQubits1QGates[i].emplace(gate.target);
        ++singleQubitMultiplicities[i][gate.target];
      } else {
        activeQubits[i].emplace(gate.control);
        activeQubits[i].emplace(gate.target);
        activeQubits2QGates[i].emplace(gate.control);
        activeQubits2QGates[i].emplace(gate.target);
        if (gate.control >= gate.target) {
          const auto edge =
              std::pair(gate.target, static_cast<std::uint16_t>(gate.control));
          if (twoQubitMultiplicities[i].find(edge) ==
              twoQubitMultiplicities[i].end()) {
            twoQubitMultiplicities[i][edge] = {0, 1};
          } else {
            twoQubitMultiplicities[i][edge].second++;
          }
        } else {
          const auto edge =
              std::pair(static_cast<std::uint16_t>(gate.control), gate.target);
          if (twoQubitMultiplicities[i].find(edge) ==
              twoQubitMultiplicities[i].end()) {
            twoQubitMultiplicities[i][edge] = {1, 0};
          } else {
            twoQubitMultiplicities[i][edge].first++;
          }
        }
      }
    }
  }
}

bool Mapper::isLayerSplittable(std::size_t index) {
  if (twoQubitMultiplicities.at(index).size() > 1) {
    return true;
  }
  if (activeQubits1QGates.at(index).size() > 2) {
    return true;
  }
  if (twoQubitMultiplicities.at(index).empty()) {
    return false;
  }
  // check if there is a 1Q gate on a qubit that is not part of the 2Q gate
  return std::any_of(activeQubits1QGates.at(index).begin(),
                     activeQubits1QGates.at(index).end(),
                     [this, index](auto q) {
                       return activeQubits2QGates.at(index).find(q) ==
                              activeQubits2QGates.at(index).end();
                     });
}

void Mapper::splitLayer(std::size_t index, Architecture& arch) {
  const SingleQubitMultiplicity& singleQubitMultiplicity =
      singleQubitMultiplicities.at(index);
  const TwoQubitMultiplicity& twoQubitMultiplicity =
      twoQubitMultiplicities.at(index);
  std::vector<Gate> layer0{};
  std::vector<Gate> layer1{};
  SingleQubitMultiplicity singleQubitMultiplicity0(arch.getNqubits(), 0);
  SingleQubitMultiplicity singleQubitMultiplicity1(arch.getNqubits(), 0);
  TwoQubitMultiplicity twoQubitMultiplicity0{};
  TwoQubitMultiplicity twoQubitMultiplicity1{};
  std::set<std::uint16_t> activeQubits0{};
  std::set<std::uint16_t> activeQubits1QGates0{};
  std::set<std::uint16_t> activeQubits2QGates0{};
  std::set<std::uint16_t> activeQubits1{};
  std::set<std::uint16_t> activeQubits1QGates1{};
  std::set<std::uint16_t> activeQubits2QGates1{};

  // 2Q-gates
  bool even = false;
  for (const auto& edge : twoQubitMultiplicity) {
    if (even) {
      twoQubitMultiplicity0.insert(edge);
      activeQubits0.emplace(edge.first.first);
      activeQubits0.emplace(edge.first.second);
      activeQubits2QGates0.emplace(edge.first.first);
      activeQubits2QGates0.emplace(edge.first.second);
    } else {
      twoQubitMultiplicity1.insert(edge);
      activeQubits1.emplace(edge.first.first);
      activeQubits1.emplace(edge.first.second);
      activeQubits2QGates1.emplace(edge.first.first);
      activeQubits2QGates1.emplace(edge.first.second);
    }
    even = !even;
  }

  // 1Q-gates
  even = true;
  for (std::size_t q = 0; q < singleQubitMultiplicity.size(); ++q) {
    if (singleQubitMultiplicity[q] == 0) {
      continue;
    }
    // if a qubit is also acted on by a 2Q-gate, put it on the same layer as
    // the 2Q-gate
    if (activeQubits2QGates0.find(static_cast<std::uint16_t>(q)) !=
        activeQubits2QGates0.end()) {
      singleQubitMultiplicity0[q] = singleQubitMultiplicity[q];
      activeQubits0.emplace(static_cast<std::uint16_t>(q));
      activeQubits1QGates0.emplace(static_cast<std::uint16_t>(q));
      continue;
    }
    if (activeQubits2QGates1.find(static_cast<std::uint16_t>(q)) !=
        activeQubits2QGates1.end()) {
      singleQubitMultiplicity1[q] = singleQubitMultiplicity[q];
      activeQubits1.emplace(static_cast<std::uint16_t>(q));
      activeQubits1QGates1.emplace(static_cast<std::uint16_t>(q));
      continue;
    }

    if (even) {
      singleQubitMultiplicity0[q] = singleQubitMultiplicity[q];
      activeQubits0.emplace(q);
      activeQubits1QGates0.emplace(q);
    } else {
      singleQubitMultiplicity1[q] = singleQubitMultiplicity[q];
      activeQubits1.emplace(q);
      activeQubits1QGates1.emplace(q);
    }
    even = !even;
  }

  for (auto& gate : layers[index]) {
    if (gate.singleQubit()) {
      if (singleQubitMultiplicity0[gate.target] > 0) {
        layer0.emplace_back(gate);
      } else {
        layer1.emplace_back(gate);
      }
    } else {
      if (activeQubits2QGates0.find(gate.target) !=
          activeQubits2QGates0.end()) {
        layer0.emplace_back(gate);
      } else {
        layer1.emplace_back(gate);
      }
    }
  }

  // insert new layers
  layers[index] = layer0;
  layers.insert(
      layers.begin() +
          static_cast<std::vector<std::vector<Mapper::Gate>>::difference_type>(
              index) +
          1,
      layer1);
  singleQubitMultiplicities[index] = singleQubitMultiplicity0;
  singleQubitMultiplicities.insert(
      singleQubitMultiplicities.begin() +
          static_cast<std::vector<SingleQubitMultiplicity>::difference_type>(
              index) +
          1,
      singleQubitMultiplicity1);
  twoQubitMultiplicities[index] = twoQubitMultiplicity0;
  twoQubitMultiplicities.insert(
      twoQubitMultiplicities.begin() +
          static_cast<std::vector<TwoQubitMultiplicity>::difference_type>(
              index) +
          1,
      twoQubitMultiplicity1);
  activeQubits[index] = activeQubits0;
  activeQubits.insert(
      activeQubits.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits1);
  activeQubits1QGates[index] = activeQubits1QGates0;
  activeQubits1QGates.insert(
      activeQubits1QGates.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits1QGates1);
  activeQubits2QGates[index] = activeQubits2QGates0;
  activeQubits2QGates.insert(
      activeQubits2QGates.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits2QGates1);
  results.input.layers = layers.size();
}

std::size_t Mapper::getNextLayer(std::size_t idx) {
  auto next = idx + 1;
  while (next < layers.size()) {
    for (const auto& gate : layers.at(next)) {
      if (!gate.singleQubit()) {
        return next;
      }
    }
    next++;
  }
  return std::numeric_limits<std::size_t>::max();
}

void Mapper::finalizeMappedCircuit() {
  // add additional qubits if the architecture contains more qubits than the
  // circuit
  if (architecture->getNqubits() > qcMapped.getNqubits()) {
    for (auto logicalQubit = qcMapped.getNqubits();
         logicalQubit < architecture->getNqubits(); ++logicalQubit) {
      auto physicalQubit = static_cast<qc::Qubit>(logicalQubit);

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logicalQubit)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physicalQubit = 0; physicalQubit < architecture->getNqubits();
             ++(physicalQubit)) {
          if (qcMapped.initialLayout.find(physicalQubit) ==
              qcMapped.initialLayout.end()) {
            break;
          }
        }
      }

      // the added logical qubits are not used in the circuit itself, so they
      // are regarded garbage
      qcMapped.addAncillaryQubit(physicalQubit, std::nullopt);
    }
  }
  // unify quantum registers
  qcMapped.unifyQuantumRegisters();

  // append measurements according to output permutation
  if (results.config.addMeasurementsToMappedCircuit) {
    qcMapped.appendMeasurementsAccordingToOutputPermutation();
  }
}

void Mapper::placeRemainingArchitectureQubits() {
  if (qc.getNqubits() < architecture->getNqubits()) {
    for (auto logical = qc.getNqubits(); logical < architecture->getNqubits();
         ++logical) {
      auto physical = static_cast<qc::Qubit>(logical);

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logical)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physical = 0; physical < architecture->getNqubits();
             ++(physical)) {
          if (qcMapped.initialLayout.find(physical) ==
              qcMapped.initialLayout.end()) {
            break;
          }
        }
      }

      qubits.at(physical) = static_cast<std::int16_t>(logical);

      // mark architecture qubit as ancillary and garbage
      qcMapped.initialLayout[physical] = static_cast<qc::Qubit>(logical);
      qcMapped.setLogicalQubitAncillary(static_cast<qc::Qubit>(logical));
      qcMapped.setLogicalQubitGarbage(static_cast<qc::Qubit>(logical));
    }
  }
}

void Mapper::preMappingOptimizations(const Configuration& config
                                     [[maybe_unused]]) {
  if (!config.preMappingOptimizations) {
    return;
  }

  // at the moment there are no pre-mapping optimizations
}

void Mapper::postMappingOptimizations(const Configuration& config) {
  if (!config.postMappingOptimizations) {
    return;
  }

  // try to cancel adjacent CNOT gates
  qc::CircuitOptimizer::cancelCNOTs(qcMapped);
}

void Mapper::countGates(decltype(qcMapped.cbegin()) it,
                        const decltype(qcMapped.cend())& end,
                        MappingResults::CircuitInfo& info) {
  for (; it != end; ++it) {
    const auto& g = *it;
    if (g->getType() == qc::Teleportation) {
      info.gates += GATES_OF_TELEPORTATION;
      continue;
    }

    if (g->isStandardOperation()) {
      if (g->getType() == qc::SWAP) {
        auto q1 = static_cast<std::uint16_t>(g->getTargets()[0]);
        auto q2 = static_cast<std::uint16_t>(g->getTargets()[1]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity += architecture->getSwapFidelityCost(q1, q2);
        }
        if (architecture->getCouplingMap().find({q1, q2}) !=
                architecture->getCouplingMap().end() &&
            architecture->getCouplingMap().find({q2, q1}) !=
                architecture->getCouplingMap().end()) {
          // bidirectional edge
          info.gates += GATES_OF_BIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
        } else {
          // unidirectional edge
          info.gates += GATES_OF_UNIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
          info.singleQubitGates += GATES_OF_DIRECTION_REVERSE;
        }
      } else if (g->getControls().empty()) {
        ++info.singleQubitGates;
        ++info.gates;
        auto q1 = static_cast<std::uint16_t>(g->getTargets()[0]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity += architecture->getSingleQubitFidelityCost(q1);
        }
      } else {
        assert(g->getType() == qc::X);
        ++info.cnots;
        ++info.gates;
        auto q1 =
            static_cast<std::uint16_t>((*(g->getControls().begin())).qubit);
        auto q2 = static_cast<std::uint16_t>(g->getTargets()[0]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity +=
              architecture->getTwoQubitFidelityCost(q1, q2);
        }
      }
      continue;
    }

    if (const auto& cg = dynamic_cast<const qc::CompoundOperation*>(g.get());
        cg != nullptr) {
      countGates(cg->cbegin(), cg->cend(), info);
    }
  }
  info.totalFidelity = std::pow(2, -info.totalLogFidelity);
}
