//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <cstdint>
#include <deque>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace na {
/**
 * @brief Class to manage the creation of layers when traversing a quantum
 * circuit.
 * @details The class uses the qc::DAG of the circuit to create layers of gates
 * that can be executed at the same time. It can be used to create the front or
 * look ahead layer.
 */

class NeutralAtomLayer {
protected:
  using DAG = std::vector<std::deque<std::unique_ptr<qc::Operation>*>>;
  using DAGIterator = std::deque<std::unique_ptr<qc::Operation>*>::iterator;
  using DAGIterators = std::vector<DAGIterator>;

  DAG dag;
  DAGIterators iterators;
  GateList gates;
  GateList mappedSingleQubitGates;
  std::vector<GateList> candidates;

  /**
   * @brief Updates the gates for the given qubits
   * @details The function iterates over the qc::DAG and updates the gates for
   * the given qubits as far es possible.
   * @param qubitsToUpdate The qubits that have been updated
   * @param commuteWith Gates the new gates should commute with
   */
  void updateByQubits(const std::set<qc::Qubit>& qubitsToUpdate);
  /**
   * @brief Updates the candidates for the given qubits
   */
  void updateCandidatesByQubits(const std::set<qc::Qubit>& qubitsToUpdate);
  /**
   * @brief Checks the candidates and add them to the gates if possible
   * @param qubitsToUpdate The qubits that have been updated
   */
  void candidatesToGates(const std::set<qc::Qubit>& qubitsToUpdate);

  // Commutation checks
  static bool commutesWithAtQubit(const GateList& layer,
                                  const qc::Operation* opPointer,
                                  const qc::Qubit& qubit);
  static bool commuteAtQubit(const qc::Operation* opPointer1,
                             const qc::Operation* opPointer2,
                             const qc::Qubit& qubit);

public:
  // Constructor
  explicit NeutralAtomLayer(DAG graph) : dag(std::move(graph)) {
    iterators.reserve(dag.size());
    candidates.reserve(dag.size());
    for (auto& i : dag) {
      auto it = i.begin();
      iterators.emplace_back(it);
      candidates.emplace_back();
    }
  }

  /**
   * @brief Returns the current layer of gates
   * @return The current layer of gates
   */
  GateList getGates() { return gates; }
  /**
   * @brief Returns a vector of the iterator indices
   * @return A copy of the current iterator indices
   */
  std::vector<uint32_t> getIteratorOffset();
  /**
   * @brief Initializes the layer by updating all qubits starting from the
   * iterators
   * @param The iterator offset to start from
   */
  void initLayerOffset(const std::vector<uint32_t>& iteratorOffset = {});
  /**
   * @brief Removes the provided gates from the current layer and update the
   * the layer depending on the qubits of the gates.
   * @param gatesToRemove Gates to remove from the current layer
   * @param commuteWith Gates the new gates should commute with
   */
  void removeGatesAndUpdate(const GateList& gatesToRemove);
  /**
   * @brief Returns the mapped single qubit gates
   * @return The mapped single qubit gates
   */
  GateList getMappedSingleQubitGates() { return mappedSingleQubitGates; }
};

} // namespace na
