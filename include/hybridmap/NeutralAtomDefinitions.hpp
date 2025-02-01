//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"

#include <cstdint>
#include <set>
#include <utility>
#include <vector>

namespace qc {
// forward declaration
class Operation;
} // namespace qc

namespace na {
// A CoordIndex corresponds to node in the SLM grid, where an atom can be placed
// (or not).
using CoordIndex = std::uint32_t;
using CoordIndices = std::vector<CoordIndex>;
// A HwQubit corresponds to an atom in the neutral atom architecture. It can be
// used as qubit or not and occupies a certain position in the architecture.
using HwQubit = uint32_t;
using HwQubits = std::set<HwQubit>;
using HwPositions [[maybe_unused]] = std::vector<HwQubits>;
// A qc::Qubit corresponds to a qubit in the quantum circuit. It can be mapped
// to a hardware qubit.
using Qubits = std::set<qc::Qubit>;

// Swaps are between hardware qc::Qubits (one of them can be unmapped).
using Swap = std::pair<HwQubit, HwQubit>;
using Swaps = std::vector<Swap>;
using WeightedSwap = std::pair<Swap, qc::fp>;
using WeightedSwaps = std::vector<WeightedSwap>;
// The distance between two hardware qubits using SWAP gates.
using SwapDistance = int32_t;
// Moves are between coordinates (the first is occupied, the second is not).
using AtomMove = std::pair<CoordIndex, CoordIndex>;

// Used to represent operations
using GateList = std::vector<const qc::Operation*>;

} // namespace na
