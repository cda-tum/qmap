//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "limits"
#include "map"
#include "random"

namespace qc {

/**
 * @brief Class that represents the hardware qubits of a neutral atom quantum
 * @details Class that represents the hardware qubits of a neutral atom quantum
 * computer. It stores the mapping from the circuit qubits to the hardware
 * qubits and the mapping from the hardware qubits to the coordinates of the
 * neutral atoms. It also stores the swap distances between the hardware
 * qubits.
 */
class HardwareQubits {
protected:
  NeutralAtomArchitecture     arch;
  Permutation                 hwToCoordIdx;
  SymmetricMatrix             swapDistances;
  std::map<HwQubit, HwQubits> nearbyQubits;
  Permutation                 initialHwPos;

  /**
   * @brief Initializes the swap distances between the hardware qubits for the
   * trivial initial layout.
   * @details Initializes the swap distances between the hardware qubits. This
   * is only valid for the trivial initial layout.
   */
  void initTrivialSwapDistances();
  /**
   * @brief Initializes the nearby qubits for each hardware qubit.
   * @details Nearby qubits are the qubits that are closer than the interaction
   * radius. Therefore they can be swapped with a single swap operation.
   */
  void initNearbyQubits();
  /**
   * @brief Computes the nearby qubits for a single hardware qubit.
   * @details Computes the nearby qubits for a single hardware qubit. This
   * function is called by initNearbyQubits(). It uses the nearby coordinates
   * of the neutral atom architecture to compute the nearby qubits.
   * @param qubit The hardware qubit for which the nearby qubits are computed.
   */
  void computeNearbyQubits(HwQubit qubit);

  /**
   * @brief Computes the swap distance between two hardware qubits.
   * @details Computes the swap distance between two hardware qubits. This
   * function is called by getSwapDistance(). It uses a breadth-first search
   * to find the shortest path between the two qubits.
   * @param q1 The first hardware qubit.
   * @param q2 The second hardware qubit.
   */
  void computeSwapDistance(HwQubit q1, HwQubit q2);

  /**
   * @brief Resets the swap distances between the hardware qubits.
   * @details Used after each shuttling operation to reset the swap distances.
   */
  void resetSwapDistances();

public:
  // Constructors
  HardwareQubits() = delete;
  HardwareQubits(const NeutralAtomArchitecture& arch,
                 InitialCoordinateMapping       initialCoordinateMapping)
      : arch(arch), swapDistances(arch.getNqubits()) {
    switch (initialCoordinateMapping) {
    case Trivial:
      for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
        hwToCoordIdx.insert({i, i});
      }
      initTrivialSwapDistances();
      break;
    case Random:
      std::vector<CoordIndex> indices(arch.getNpositions());
      std::iota(indices.begin(), indices.end(), 0);
      std::random_device rd;
      std::mt19937       g(rd());
      std::shuffle(indices.begin(), indices.end(), g);
      for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
        hwToCoordIdx.insert({i, indices[i]});
      }

      swapDistances = SymmetricMatrix(arch.getNqubits(), -1);
    }
    initNearbyQubits();
    initialHwPos = hwToCoordIdx;
  }

  // Mapping

  /**
   * @brief Checks if a hardware qubit is mapped to a coordinate.
   * @param idx The coordinate index.
   * @return Boolean indicating if the hardware qubit is mapped to a coordinate.
   */
  [[nodiscard]] inline bool isMapped(CoordIndex idx) const {
    return !std::none_of(
        hwToCoordIdx.begin(), hwToCoordIdx.end(),
        [idx](const auto& pair) { return pair.second == idx; });
  }
  /**
   * @brief Updates mapping after moving a hardware qubit to a coordinate.
   * @details Checks if the coordinate is valid and free. If yes, the mapping is
   * updated.
   * @param hwQubit The hardware qubit to be moved.
   * @param newCoord The new coordinate of the hardware qubit.
   */
  void move(HwQubit hwQubit, CoordIndex newCoord);

  /**
   * @brief Converts gate qubits from hardware qubits to coordinate indices.
   * @param op The operation.
   */
  inline void mapToCoordIdx(Operation* op) const {
    op->setTargets(hwToCoordIdx.apply(op->getTargets()));
    op->setControls(hwToCoordIdx.apply(op->getControls()));
  }

  /**
   * @brief Returns the coordinate index of a hardware qubit.
   * @param qubit The hardware qubit.
   * @return The coordinate index of the hardware qubit.
   */
  [[nodiscard]] inline CoordIndex getCoordIndex(HwQubit qubit) const {
    return hwToCoordIdx.at(qubit);
  }
  /**
   * @brief Returns the coordinate indices of a set of hardware qubits.
   * @param hwQubits The set of hardware qubits.
   * @return The coordinate indices of the hardware qubits.
   */
  [[nodiscard]] inline std::set<CoordIndex>
  getCoordIndices(std::set<HwQubit>& hwQubits) const {
    std::set<CoordIndex> coordIndices;
    for (auto const& hwQubit : hwQubits) {
      coordIndices.insert(this->getCoordIndex(hwQubit));
    }
    return coordIndices;
  }
  /**
   * @brief Returns the hardware qubit at a coordinate.
   * @details Returns the hardware qubit at a coordinate. Throws an exception if
   * there is no hardware qubit at the coordinate.
   * @param coordIndex The coordinate index.
   * @return The hardware qubit at the coordinate.
   */
  [[nodiscard]] inline HwQubit getHwQubit(CoordIndex coordIndex) const {
    for (auto const& [hwQubit, index] : hwToCoordIdx) {
      if (index == coordIndex) {
        return hwQubit;
      }
    }
    throw std::runtime_error("There is no qubit at this coordinate " +
                             std::to_string(coordIndex));
  }

  // Forwards from architecture class

  /**
   * @brief Returns the nearby coordinates of a hardware qubit.
   * @param q The hardware qubit.
   * @return The nearby coordinates of the hardware qubit.
   */
  [[nodiscard]] inline std::set<CoordIndex>
  getNearbyCoordinates(HwQubit q) const {
    return this->arch.getNearbyCoordinates(this->getCoordIndex(q));
  }

  // Swap Distances and Nearby Qubits

  /**
   * @brief Returns the swap distance between two hardware qubits.
   * @details Returns the swap distance between two hardware qubits. If the
   * swap distance is not yet computed, it is computed using a breadth-first
   * search.
   * @param q1 The first hardware qubit.
   * @param q2 The second hardware qubit.
   * @param closeBy If the swap should be performed to the exact position of q2
   * or just to its vicinity.
   * @return The swap distance between the two hardware qubits.
   */
  [[nodiscard]] inline fp getSwapDistance(HwQubit q1, HwQubit q2,
                                          bool closeBy = true) {
    if (q1 == q2) {
      return 0;
    }
    if (swapDistances(q1, q2) < 0) {
      computeSwapDistance(q1, q2);
    }
    if (closeBy) {
      return swapDistances(q1, q2);
    }
    return swapDistances(q1, q2) + 1;
  }
  /**
   * @brief Returns the swap distance from a hardware qubit to a coordinate.
   * @details Checks the vicinity of the coordinate for a free place and returns
   * the minimal swap distance of all free places. Returns infinity if there is
   * no free place in the vicinity.
   * @param idx The coordinate index.
   * @param target The hardware qubit.
   * @return The swap distance from the hardware qubit to the coordinate.
   */
  [[nodiscard]] fp getSwapDistanceMove(CoordIndex idx, HwQubit target);

  /**
   * @brief Returns the nearby hardware qubits of a hardware qubit.
   * @param q The hardware qubit.
   * @return The nearby hardware qubits of the hardware qubit.
   */
  [[nodiscard]] inline HwQubits getNearbyQubits(HwQubit q) const {
    return nearbyQubits.at(q);
  }

  /**
   * @brief Returns vector of all possible swaps for a hardware qubit.
   * @param q The hardware qubit.
   * @return The vector of all possible swaps for the hardware qubit.
   */
  std::vector<Swap> getNearbySwaps(HwQubit q) const;
  /**
   * @brief Returns the unoccupied coordinates in the vicinity of a hardware
   * qubit.
   * @param q The hardware qubit.
   * @return The unoccupied coordinates in the vicinity of the hardware qubit.
   */
  std::set<CoordIndex> getNearbyFreeCoordinates(HwQubit q);

  /**
   * @brief Returns the unoccupied coordinates in the vicinity of a coordinate.
   * @param idx The coordinate index.
   * @return The unoccupied coordinates in the vicinity of the coordinate.
   */
  std::set<CoordIndex> getNearbyFreeCoordinatesByCoord(CoordIndex idx);
  /**
   * @brief Returns the occupied coordinates in the vicinity of a hardware
   * qubit.
   * @param q The hardware qubit.
   * @return The occupied coordinates in the vicinity of the hardware qubit.
   */
  std::set<CoordIndex> getNearbyOccupiedCoordinates(HwQubit q);

  /**
   * @brief Returns the occupied coordinates in the vicinity of a coordinate.
   * @param idx The coordinate index.
   * @return The occupied coordinates in the vicinity of the coordinate.
   */
  std::set<CoordIndex> getNearbyOccupiedCoordinatesByCoord(CoordIndex idx);

  /**
   * @brief Computes the summed swap distance between all hardware qubits in a
   * set.
   * @param qubits The set of hardware qubits.
   * @return The summed swap distance between all hardware qubits in the set.
   */
  fp getAllToAllSwapDistance(std::set<HwQubit>& qubits);

  /**
   * @brief Computes the closest free coordinate in a given direction.
   * @details Uses a breadth-first search to find the closest free coordinate in
   * a given direction.
   * @param qubit The hardware qubit to start the search from.
   * @param direction The direction in which the search is performed
   * (Left/Right, Down/Up)
   * @param excludedCoords Coordinates to be ignored in the search.
   * @return The closest free coordinate in the given direction.
   */
  std::vector<CoordIndex>
  findClosestFreeCoord(HwQubit qubit, Direction direction,
                       const CoordIndices& excludedCoords = {});

  // Blocking
  /**
   * @brief Computes all hardware qubits that are blocked by a set of hardware
   * qubits.
   * @param qubits The input hardware qubits.
   * @return The blocked hardware qubits.
   */
  std::set<HwQubit> getBlockedQubits(const std::set<HwQubit>& qubits);

  std::map<HwQubit, HwQubit> getInitialHwPos() const {
    std::map<HwQubit, HwQubit> initialHwPosMap;
    for (auto const& pair : initialHwPos) {
      initialHwPosMap[pair.first] = pair.second;
    }
    return initialHwPosMap;
  }
};
} // namespace qc