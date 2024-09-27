//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/Tableau.hpp"

#include "ir/QuantumComputation.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <istream>
#include <optional>
#include <ostream>
#include <plog/Log.h>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cs {
void Tableau::dump(const std::string& filename) const {
  auto of = std::ofstream(filename);
  if (!of.good()) {
    const auto msg = "Error opening file " + filename;
    PLOG_FATAL << msg;
    throw std::runtime_error(msg);
  }
  dump(of);
}

void Tableau::dump(std::ostream& of) const { of << *this; }

void Tableau::import(const std::string& filename) {
  auto is = std::ifstream(filename);
  if (!is.good()) {
    const auto msg = "Error opening file " + filename;
    PLOG_FATAL << msg;
    throw std::runtime_error(msg);
  }
  import(is);
}

void parseLine(const std::string& line, char separator,
               const std::set<char>& escapeChars,
               const std::set<char>& ignoredChars,
               std::vector<std::string>& result) {
  result.clear();
  std::string word;
  bool inEscape = false;
  for (const char c : line) {
    if (ignoredChars.find(c) != ignoredChars.end()) {
      continue;
    }
    if (inEscape) {
      if (escapeChars.find(c) != escapeChars.end()) {
        inEscape = false;
      } else {
        word += c;
      }
    } else {
      if (escapeChars.find(c) != escapeChars.end()) {
        inEscape = true;
      } else if (c == separator) {
        result.push_back(word);
        word = "";
      } else {
        word += c;
      }
    }
  }
  result.push_back(word);
}

void Tableau::import(std::istream& is) {
  tableau.clear();

  std::string line;
  std::vector<std::string> data{};
  char delimiter = '|';

  while (std::getline(is, line)) {
    if (line.find('|', 0) == std::string::npos) {
      delimiter = ';';
    }
    tableau.emplace_back();
    parseLine(line, delimiter, {'\"'}, {'\\', '\r', '\n', '\t'}, data);
    for (const auto& datum : data) {
      if (datum.empty()) {
        continue;
      }
      tableau.back().emplace_back(static_cast<EntryType>(std::stoul(datum)));
    }
  }
}

void Tableau::applyGate(const qc::Operation* const gate) {
  if (gate->getNcontrols() > 1U) {
    const auto* const msg =
        "Tableau::applyGate: Only operations with up to one control "
        "are supported.";
    PLOG_FATAL << msg;
    throw std::runtime_error(msg);
  }

  const auto target = static_cast<std::size_t>(gate->getTargets().at(0U));

  // non-controlled gates
  if (!gate->isControlled()) {
    switch (gate->getType()) {
    case qc::OpType::H:
      applyH(target);
      break;
    case qc::OpType::S:
      applyS(target);
      break;
    case qc::OpType::Sdg:
      applySdag(target);
      break;
    case qc::OpType::SX:
      applySx(target);
      break;
    case qc::OpType::SXdg:
      applySxdag(target);
      break;
    case qc::OpType::X:
      applyX(target);
      break;
    case qc::OpType::Y:
      applyY(target);
      break;
    case qc::OpType::Z:
      applyZ(target);
      break;
    case qc::OpType::SWAP: {
      const auto target2 = static_cast<std::size_t>(gate->getTargets().at(1U));
      applySwap(target, target2);
      break;
    }
    case qc::OpType::iSWAP: {
      const auto target2 = static_cast<std::size_t>(gate->getTargets().at(1U));
      applyISwap(target, target2);
      break;
    }
    case qc::OpType::DCX: {
      const auto target2 = static_cast<std::size_t>(gate->getTargets().at(1U));
      applyDCX(target, target2);
      break;
    }
    case qc::OpType::ECR: {
      const auto target2 = static_cast<std::size_t>(gate->getTargets().at(1U));
      applyECR(target, target2);
      break;
    }
    default:
      // unsupported non-controlled gate type
      const auto msg =
          "Tableau::applyGate: Unsupported non-controlled gate type " +
          qc::toString(gate->getType());
      PLOG_FATAL << msg;
      throw std::runtime_error(msg);
    }
  } else {
    const auto control =
        static_cast<std::size_t>((*gate->getControls().begin()).qubit);
    switch (gate->getType()) {
    case qc::OpType::X:
      applyCX(control, target);
      break;
    case qc::OpType::Y:
      applyCY(control, target);
      break;
    case qc::OpType::Z:
      applyCZ(control, target);
      break;
    default:
      // unsupported controlled gate type
      const auto msg = "Tableau::applyGate: Unsupported controlled gate type " +
                       qc::toString(gate->getType());
      PLOG_FATAL << msg;
      throw std::runtime_error(msg);
    }
  }
}

void Tableau::createDiagonalTableau(const std::size_t nQ,
                                    const bool includeDestabilizers) {
  nQubits = nQ;
  tableau.clear();
  if (includeDestabilizers) {
    tableau.resize(2U * nQubits);
  } else {
    tableau.resize(nQubits);
  }
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i].resize((2U * nQubits) + 1U);
    if (includeDestabilizers) {
      for (std::size_t j = 0; j < (2U * nQubits); ++j) {
        tableau[i][j] = 0;
        if (i == j) {
          tableau[i][j] = 1;
        }
      }
    } else {
      for (std::size_t j = nQubits; j < (2U * nQubits); ++j) {
        tableau[i][j] = 0;
        if (i == (j - nQubits)) {
          tableau[i][j] = 1;
        }
      }
    }
  }
}

std::string Tableau::toString() const {
  std::stringstream ss;
  for (const auto& row : tableau) {
    if (row.size() != tableau.back().size()) {
      const auto* const msg = "Tableau::toString: Tableau is not rectangular";
      PLOG_FATAL << msg;
      throw std::runtime_error(msg);
    }
    for (const auto& s : row) {
      ss << std::to_string(s) << ';';
    }
    ss << "\n";
  }
  return ss.str();
}

void Tableau::fromString(const std::string& str) {
  std::stringstream ss(str);
  std::string line;
  std::getline(ss, line);
  if (line.empty()) {
    return;
  }
  const auto rStabilizer = std::regex("([\\+-]?)([IYZX]+)");
  if (std::regex_search(line, rStabilizer)) {
    // string is a list of stabilizers
    loadStabilizerDestabilizerString(str);
  } else {
    // assume string is a semicolon separated binary matrix
    ss = std::stringstream(str);
    import(ss);
  }
}

void Tableau::applyH(const std::size_t target) {
  assert(target < nQubits);
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i][2U * nQubits] ^= static_cast<EntryType>(
        tableau[i][target] & tableau[i][target + nQubits]);
    std::swap(tableau[i][target], tableau[i][target + nQubits]);
  }
}

void Tableau::applyS(const std::size_t target) {
  assert(target < nQubits);
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i][2U * nQubits] ^= static_cast<EntryType>(
        tableau[i][target] & tableau[i][target + nQubits]);
    tableau[i][target + nQubits] ^= tableau[i][target];
  }
}

// Sdag = S * S * S
void Tableau::applySdag(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyS(target);
  applyS(target);
}

// Sx = Sdag * H * Sdag
void Tableau::applySx(const std::size_t target) {
  assert(target < nQubits);
  applySdag(target);
  applyH(target);
  applySdag(target);
}

// Sxdag = S * H * S
void Tableau::applySxdag(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyH(target);
  applyS(target);
}

// X = H * Z * H
void Tableau::applyX(const std::size_t target) {
  assert(target < nQubits);
  applyH(target);
  applyZ(target);
  applyH(target);
}

// Y = X * Z
void Tableau::applyY(const std::size_t target) {
  assert(target < nQubits);
  applyX(target);
  applyZ(target);
}

// Z = S * S
void Tableau::applyZ(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyS(target);
}

void Tableau::applyCX(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  for (auto i = 0U; i < getTableauSize(); ++i) {
    const auto xa = tableau[i][control];
    const auto za = tableau[i][control + nQubits];
    const auto xb = tableau[i][target];
    const auto zb = tableau[i][target + nQubits];
    tableau[i][2 * nQubits] ^= (xa & zb) & ((xb ^ za) ^ 1);
    tableau[i][control + nQubits] = za ^ zb;
    tableau[i][target] = xb ^ xa;
  }
}

void Tableau::applyCY(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  applySdag(target);
  applyCX(control, target);
  applyS(target);
}

void Tableau::applyCZ(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  applyH(target);
  applyCX(control, target);
  applyH(target);
}

void Tableau::applySwap(const std::size_t q1, const std::size_t q2) {
  assert(q1 < nQubits);
  assert(q2 < nQubits);
  assert(q1 != q2);
  applyCX(q1, q2);
  applyCX(q2, q1);
  applyCX(q1, q2);
}

void Tableau::applyISwap(const std::size_t q1, const std::size_t q2) {
  assert(q1 < nQubits);
  assert(q2 < nQubits);
  assert(q1 != q2);
  applyS(q2);
  applyS(q1);
  applyH(q1);
  applyDCX(q1, q2);
  applyH(q2);
}

void Tableau::applyDCX(const std::size_t q1, const std::size_t q2) {
  assert(q1 < nQubits);
  assert(q2 < nQubits);
  assert(q1 != q2);
  applyCX(q1, q2);
  applyCX(q2, q1);
}

void Tableau::applyECR(const std::size_t q1, const std::size_t q2) {
  assert(q1 < nQubits);
  assert(q2 < nQubits);
  assert(q1 != q2);
  applyS(q1);
  applySx(q2);
  applyCX(q1, q2);
  applyX(q1);
}

Tableau::Tableau(const qc::QuantumComputation& qc, const std::size_t begin,
                 const std::size_t end, const bool includeDestabilizers)
    : Tableau(qc.getNqubits(), includeDestabilizers) {
  std::size_t currentG = 0;
  for (const auto& gate : qc) {
    if (gate->getType() == qc::OpType::Compound) {
      const auto* const compOp =
          dynamic_cast<const qc::CompoundOperation* const>(gate.get());
      auto cit = compOp->begin();
      while (cit != compOp->end()) {
        if (currentG >= begin && (currentG < end)) {
          applyGate((*cit).get());
        }
        ++cit;
        ++currentG;
      }
    } else {
      if (currentG >= begin && (currentG < end)) {
        applyGate(gate.get());
      }
      ++currentG;
    }
    if (currentG >= end) {
      break;
    }
  }
}
void Tableau::fromString(const std::string& stabilizers,
                         const std::string& destabilizers) {
  loadStabilizerDestabilizerString(destabilizers);
  loadStabilizerDestabilizerString(stabilizers);
}

Tableau::RowType Tableau::parseStabilizer(const std::string& stab) {
  auto stabCopy = stab;

  if (stabCopy[0] == '\'') {
    if (stabCopy[stabCopy.size() - 1] == '\'') {
      stabCopy = stabCopy.substr(1, stabCopy.size() - 2);
    } else {
      throw std::invalid_argument("Unmatched \"'\" in stabilizer string");
    }
  }
  if (stabCopy[0] == '+' || stabCopy[0] == '-') {
    stabCopy = stabCopy.substr(1);
  }

  RowType row;
  for (const auto c : stabCopy) {
    if (c == 'I' || c == '_' || c == 'Z') {
      row.push_back(0);
    } else if (c == 'X' || c == 'Y') {
      row.push_back(1);
    } else {
      throw std::invalid_argument(
          "Invalid stabilizer " + stab +
          ". Stabilizers must be given as a list of stabilizer "
          "like [XYZI, ZIXZ]");
    }
  }
  for (const auto c : stabCopy) {
    if (c == 'I' || c == '_' || c == 'X') {
      row.push_back(0);
    } else if (c == 'Y' || c == 'Z') {
      row.push_back(1);
    } else {
      throw std::invalid_argument(
          "Invalid stabilizer" + stab +
          ". Stabilizers must be given as a list of stabilizer "
          "like [XYZI, ZIXZ]");
    }
  }
  if (stab[0U] == '-' || (stab[0U] == '\'' && stab[1U] == '-')) {
    row.push_back(1);
  } else {
    row.push_back(0);
  }
  return row;
}

void Tableau::loadStabilizerDestabilizerString(const std::string& string) {
  std::stringstream ss(string);
  std::string line;
  std::getline(ss, line);
  if (line.empty()) {
    return;
  }

  auto stabilizers = line;
  stabilizers.erase(remove_if(stabilizers.begin(), stabilizers.end(), isspace),
                    stabilizers.end());

  if (stabilizers[0] == '[') {
    if (stabilizers[stabilizers.size() - 1] == ']') {
      stabilizers = stabilizers.substr(1, stabilizers.size() - 2);
    } else {
      throw std::invalid_argument("Unmatched \"[\" in stabilizer string");
    }
  }

  std::optional<std::size_t> stabLength;
  const auto& checkStabLength = [&](const RowType& row) {
    if (!stabLength.has_value()) {
      stabLength = row.size();
    }
    if (stabLength.value() != row.size()) {
      throw std::invalid_argument("All Stabilizers muts have the same length");
    }
  };

  const char delimiter = ',';
  std::string stab;
  for (std::size_t pos = stabilizers.find(delimiter); pos != std::string::npos;
       pos = stabilizers.find(delimiter)) {
    stab = stabilizers.substr(0, pos);
    const auto& row = parseStabilizer(stab);
    checkStabLength(row);
    tableau.push_back(row);
    stabilizers = stabilizers.substr(pos + 1);
  }
  const auto& row =
      parseStabilizer(stabilizers); // parse stabilizer past last comma
  checkStabLength(row);
  tableau.push_back(row);
}
bool Tableau::isIdentityTableau() const {
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    for (std::size_t j = 0U; j < tableau[i].size(); ++j) {
      if ((j == i && tableau[i][j] == 0U) || (j != i && tableau[i][j] == 1U)) {
        return false;
      }
    }
  }
  return true;
}

Tableau Tableau::applyMapping(const std::vector<std::vector<bool>>& p) {
  Tableau mappedTableau = Tableau(nQubits, hasDestabilizers());
  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < p[i].size(); j++) {
      // apply mapping from column i to j if p is set
      if (p[i][j]) {
        // in every row swap x entry and z entry
        for (size_t n = 0; n < mappedTableau.getTableauSize(); n++) {
          mappedTableau.tableau[n][j] = tableau[n][i];
          mappedTableau.tableau[n][j + mappedTableau.nQubits] =
              tableau[n][i + mappedTableau.nQubits];
        }
      }
    }
  }
  // copy r column without changes
  for (size_t i = 0; i < tableau.size(); i++) {
    mappedTableau.tableau[i][2 * nQubits] = tableau[i][2 * nQubits];
  }
  return mappedTableau;
}

// number of Qubits is passed because nQubits is not set in result Tableau of
// synthesis
Tableau Tableau::reverseMappingOnRows(const std::vector<std::vector<bool>>& p,
                                      const size_t nq) {
  auto mappedTableau = Tableau(nq, true);
  mappedTableau.tableau = tableau;
  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < p[i].size(); j++) {
      // apply mapping from row i to j if p is set
      if (p[i][j]) {
        mappedTableau.tableau[i] = tableau[j];
        mappedTableau.tableau[i + nq] = tableau[j + nq];
      }
    }
  }
  return mappedTableau;
}

// in place Gauss Elimination
void Tableau::rref() {
  const size_t rows = tableau.size();
  const size_t cols = tableau[0].size();
  if (rows == 1) {
    return;
  }

  size_t pivot = 0;

  for (size_t col = 0; col < cols; ++col) {
    // find the pivot row for the current column
    size_t maxRow = pivot;
    if (maxRow >= rows) {
      break;
    }
    for (size_t row = pivot; row < rows; ++row) {
      if (tableau[row][col] == 1) {
        maxRow = row;
        break;
      }
    }
    // if no pivot is found, continue to the next column
    if (tableau[maxRow][col] == 0) {
      continue;
    }
    // swap the pivot row with the current row
    std::swap(tableau[pivot], tableau[maxRow]);

    // eliminate all other 1s in the current column
    for (size_t row = 0; row < rows; ++row) {
      if (row != pivot && tableau[row][col] == 1) {
        for (size_t k = 0; k < cols; ++k) {
          if (tableau[row][k] == tableau[pivot][k]) {
            tableau[row][k] = 0;
          } else {
            tableau[row][k] = 1;
          }
        }
      }
    }
    pivot++;
  }
}
// equivalence check for two Tableaus without destablilizers
bool Tableau::equivalentUpToStabilizer(const Tableau* t) const {
  if (getTableauSize() != t->getTableauSize()) {
    return false;
  }

  Tableau tEliminated = Tableau(nQubits, hasDestabilizers());
  tEliminated.tableau = tableau;
  tEliminated.rref();

  Tableau t2Eliminated = Tableau(nQubits, hasDestabilizers());
  t2Eliminated.tableau = t->getTableau();
  t2Eliminated.rref();
  return tEliminated == t2Eliminated;
}
} // namespace cs
