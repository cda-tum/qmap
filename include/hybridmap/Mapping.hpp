//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"

namespace qc {

/**
 * @brief Class to manage the mapping between circuit qubits and hardware qubits
 * in a bijective manner.
 */
class Mapping {
protected:
  // std::map<Qubit, HwQubit>
  Permutation circToHw;

public:
  Mapping() = default;
  Mapping(size_t nQubits, InitialMapping initialMapping) {
    switch (initialMapping) {
    case Identity:
      for (size_t i = 0; i < nQubits; ++i) {
        circToHw.insert({i, i});
      }
      break;
    }
  }
  /**
   * @brief Assigns a circuit qubit to a hardware qubit.
   * @param qubit The circuit qubit to be assigned
   * @param hwQubit The hardware qubit to be assigned
   */
  void inline setCircuitQubit(Qubit qubit, HwQubit hwQubit) {
    circToHw[qubit] = hwQubit;
  }

  /**
   * @brief Returns the hardware qubit assigned to the given circuit qubit.
   * @param qubit The circuit qubit to be queried
   * @return The hardware qubit assigned to the given circuit qubit
   */
  [[nodiscard]] inline HwQubit getHwQubit(Qubit qubit) const {
    return circToHw.at(qubit);
  }

  /**
   * @brief Returns the hardware qubits assigned to the given circuit qubits.
   * @param qubits The circuit qubits to be queried
   * @return The hardware qubits assigned to the given circuit qubits
   */
  [[nodiscard]] inline std::set<HwQubit>
  getHwQubits(std::set<Qubit>& qubits) const {
    std::set<HwQubit> hwQubits;
    for (const auto& qubit : qubits) {
      hwQubits.insert(this->getHwQubit(qubit));
    }
    return hwQubits;
  }

  /**
   * @brief Returns the circuit qubit assigned to the given hardware qubit.
   * @details Throws an exception if the hardware qubit is not assigned to any
   * circuit qubit.
   * @param qubit The hardware qubit to be queried
   * @return The circuit qubit assigned to the given hardware qubit
   */
  [[nodiscard]] inline Qubit getCircQubit(HwQubit qubit) const {
    for (const auto& [circQubit, hwQubit] : circToHw) {
      if (hwQubit == qubit) {
        return circQubit;
      }
    }
    throw std::runtime_error("Hardware qubit: " + std::to_string(qubit) +
                             " not found in mapping");
  }

  /**
   * @brief Indicates if any circuit qubit is assigned to the given hardware
   * qubit.
   * @param qubit The hardware qubit to be queried
   * @return True if any circuit qubit is assigned to the given hardware qubit,
   * false otherwise
   */
  [[nodiscard]] inline bool isMapped(HwQubit qubit) const {
    return std::any_of(
        circToHw.begin(), circToHw.end(),
        [qubit](const auto& pair) { return pair.second == qubit; });
  }

  /**
   * @brief Converts the qubits of an operation from circuit qubits to hardware
   * qubits.
   * @param op The operation to be converted
   */
  inline void mapToHwQubits(Operation* op) const {
    op->setTargets(circToHw.apply(op->getTargets()));
    op->setControls(circToHw.apply(op->getControls()));
  }

  /**
   * @brief Interchanges the mapping of two hardware qubits. At least one of it
   * must be mapped to a circuit qubit.
   * @param swap The two circuit qubits to be swapped
   */
  void swap(Swap swap);
};

} // namespace qc