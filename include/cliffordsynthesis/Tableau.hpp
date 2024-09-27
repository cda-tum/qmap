//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"

#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cs {
class Tableau {
  using EntryType = std::uint8_t;
  using RowType = std::vector<EntryType>;
  using TableauType = std::vector<RowType>;
  std::size_t nQubits{};
  TableauType tableau;

private:
  void loadStabilizerDestabilizerString(const std::string& string);
  static RowType parseStabilizer(const std::string& stab);

public:
  Tableau() = default;
  explicit Tableau(const qc::QuantumComputation& qc, std::size_t begin = 0,
                   std::size_t end = std::numeric_limits<std::size_t>::max(),
                   bool includeDestabilizers = false);
  explicit Tableau(const std::size_t nq,
                   const bool includeDestabilizers = false)
      : nQubits(nq) {
    createDiagonalTableau(nq, includeDestabilizers);
  }
  explicit Tableau(const std::string& description) {
    fromString(description);
    if (tableau.empty()) {
      throw std::runtime_error("Tableau is empty");
    }
    nQubits = tableau.back().size() / 2U;
  }
  explicit Tableau(const std::string& stabilizers,
                   const std::string& destabilizers) {
    fromString(stabilizers, destabilizers);
    nQubits = tableau.size() / 2U;
  }

  [[nodiscard]] RowType operator[](const std::size_t index) {
    return tableau[index];
  }
  [[nodiscard]] RowType operator[](const std::size_t index) const {
    return tableau[index];
  }

  [[nodiscard]] RowType at(const std::size_t index) {
    return tableau.at(index);
  }

  [[nodiscard]] std::size_t getQubitCount() const { return nQubits; }

  [[nodiscard]] std::size_t getTableauSize() const { return tableau.size(); }

  [[nodiscard]] bool hasDestabilizers() const {
    return tableau.size() == 2 * nQubits;
  }

  [[nodiscard]] auto& getTableau() const { return tableau; }

  void dump(const std::string& filename) const;

  void dump(std::ostream& of) const;

  void import(const std::string& filename);
  void import(std::istream& is);

  template <std::size_t N>
  void populateTableauFrom(const std::bitset<N> bv, const std::size_t nQ,
                           const std::size_t column) {
    assert(column <= 2 * nQ);
    assert(nQ <= getTableauSize());
    assert(nQ <= N);
    for (std::size_t i = 0U; i < nQ; ++i) {
      if (bv[i]) {
        tableau[i][column] = 1;
      } else {
        tableau[i][column] = 0;
      }
    }
  }
  void populateTableauFrom(const std::uint64_t bv, const std::size_t nQ,
                           const std::size_t column) {
    populateTableauFrom(std::bitset<64>(bv), nQ, column);
  }

  void applyGate(const qc::Operation* gate);
  void applyH(std::size_t target);
  void applyS(std::size_t target);
  void applySdag(std::size_t target);
  void applySx(std::size_t target);
  void applySxdag(std::size_t target);
  void applyX(std::size_t target);
  void applyY(std::size_t target);
  void applyZ(std::size_t target);
  void applyCX(std::size_t control, std::size_t target);
  void applyCY(std::size_t control, std::size_t target);
  void applyCZ(std::size_t control, std::size_t target);
  void applySwap(std::size_t q1, std::size_t q2);
  void applyISwap(std::size_t q1, std::size_t q2);
  void applyDCX(std::size_t q1, std::size_t q2);
  void applyECR(std::size_t q1, std::size_t q2);

  [[gnu::pure]] friend bool operator==(const Tableau& lhs, const Tableau& rhs) {
    return lhs.tableau == rhs.tableau;
  }
  [[gnu::pure]] friend bool operator!=(const Tableau& lhs, const Tableau& rhs) {
    return !(lhs == rhs);
  }

  [[nodiscard]] bool isIdentityTableau() const;

  void createDiagonalTableau(std::size_t nQ, bool includeDestabilizers = false);

  friend std::ostream& operator<<(std::ostream& os, const Tableau& dt) {
    os << dt.toString();
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Tableau& dt) {
    if (is.good()) {
      std::stringstream ss;
      ss << is.rdbuf();
      dt.fromString(ss.str());
    }
    return is;
  }

  [[nodiscard]] std::string toString() const;
  void fromString(const std::string& str);

  void fromString(const std::string& stabilizers,
                  const std::string& destabilizers);
  template <std::size_t N>
  [[nodiscard]] std::bitset<N> getBVFrom(const std::size_t column) const {
    assert(column <= 2 * nQubits);
    assert(nQubits <= N);
    std::bitset<N> bv;
    for (std::size_t i = 0U; i < getTableauSize(); ++i) {
      if (tableau[i][column] == 1U) {
        bv[i] = 1;
      }
    }
    return bv;
  }
  [[nodiscard]] std::uint64_t getBVFrom(const std::size_t column) const {
    return getBVFrom<64>(column).to_ullong();
  }
  Tableau applyMapping(const std::vector<std::vector<bool>>& p);
  Tableau reverseMappingOnRows(const std::vector<std::vector<bool>>& p,
                               size_t nq);
  void rref();
  bool equivalentUpToStabilizer(const Tableau* t) const;
};
} // namespace cs
