#include "LogicTerm.hpp"

#include "Logic.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace logicbase {

uint64_t LogicTerm::getMaxChildrenDepth() const {
  uint64_t max = 0;
  for (const LogicTerm& t : getNodes()) {
    const uint64_t d = t.getMaxChildrenDepth();
    if (d > max) {
      max = d;
    }
  }
  return max + 1;
}

std::string LogicTerm::getStrRep(OpType op) {
  std::stringstream os;
  switch (op) {
  case OpType::Constant:
    os << "CONST";
    break;
  case OpType::Variable:
    os << "VAR";
    break;
  case OpType::AND:
    os << "<AND ";
    break;
  case OpType::OR:
    os << "<OR ";
    break;
  case OpType::BitAnd:
    os << "<BV_AND ";
    break;
  case OpType::BitOr:
    os << "<BV_OR ";
    break;
  case OpType::ITE:
    os << "<ITE ";
    break;
  case OpType::NEG:
    os << "<NEG ";
    break;
  case OpType::EQ:
    os << "<EQ ";
    break;
  case OpType::XOR:
    os << "<XOR ";
    break;
  case OpType::BitEq:
    os << "<BV_EQ ";
    break;
  case OpType::BitXor:
    os << "<BV_XOR ";
    break;
  case OpType::IMPL:
    os << "<IMPL ";
    break;
  case OpType::ADD:
    os << "<ADD ";
    break;
  case OpType::SUB:
    os << "<SUB ";
    break;
  case OpType::MUL:
    os << "<MUL ";
    break;
  case OpType::DIV:
    os << "<DIV ";
    break;
  case OpType::GT:
    os << "<GT ";
    break;
  case OpType::LT:
    os << "<LT ";
    break;
  case OpType::GTE:
    os << "<GTE ";
    break;
  case OpType::LTE:
    os << "<LTE ";
    break;
  default:
    os << "<ERROR TYPE";
    break;
  }
  return os.str();
}

LogicTerm::LogicTerm(const OpType op, const LogicTerm& a, const LogicTerm& b)
    : lb(getValidLogicPtr(a, b)) {
  if (a.isConst() || b.isConst()) {
    *this = combineConst(a, b, op, lb);
    return;
  }
  if (op == OpType::AND || op == OpType::OR) {
    *this = combineTerms(a.getBoolConversion(), b.getBoolConversion(), op, lb);
    return;
  }
  *this = combineTerms(a, b, op, lb);
}

LogicTerm::LogicTerm(OpType op, const LogicTerm& a) {
  if (op == OpType::NEG) {
    *this = negateTerm(a, a.getLogic());
    return;
  }
  throw std::runtime_error("Invalid opType");
}

LogicTerm::LogicTerm(OpType op, const LogicTerm& a, const LogicTerm& b,
                     const LogicTerm& c)
    : LogicTerm(op, a, b, c, getTargetCType(b, c, op),
                getValidLogicPtr(a, b, c)) {}

uint64_t LogicTerm::getNextId(Logic* logic) {
  if (logic == nullptr) {
    return gid++;
  }
  return logic->getNextId();
}

LogicTerm::LogicTerm(OpType op, const std::initializer_list<LogicTerm>& n,
                     CType type, Logic* logic)
    : lb(logic), id(getNextId(logic)), depth(getMax(n)), name(getStrRep(op)),
      opType(op), bvSize(getMaxBVSize(n)), nodes(n), cType(type) {}

LogicTerm::LogicTerm(OpType op, const std::vector<LogicTerm>& n, CType type,
                     Logic* logic)
    : lb(logic), id(getNextId(logic)), depth(getMax(n)), name(getStrRep(op)),
      opType(op), bvSize(getMaxBVSize(n)), nodes(n), cType(type) {}

LogicTerm::LogicTerm(Logic* logic)
    : lb(logic), id(getNextId(lb)), name(std::to_string(id)) {}

LogicTerm::LogicTerm(std::string n, Logic* logic)
    : lb(logic), id(getNextId(lb)), name(std::move(n)) {}

LogicTerm::LogicTerm(OpType op, std::string n, CType type,
                     Logic* logic) // potentially , uint16_t bvs = 0
    : lb(logic), id(getNextId(lb)), name(std::move(n)), opType(op),
      cType(type) {}

LogicTerm::LogicTerm(std::string n, const uint64_t identifier, Logic* logic)
    : lb(logic), id(identifier), name(std::move(n)) {}

LogicTerm::LogicTerm(CType type, Logic* logic)
    : lb(logic), id(getNextId(lb)), name(std::to_string(id)), cType(type) {}

LogicTerm::LogicTerm(std::string n, CType type, Logic* logic, uint16_t bvs)
    : lb(logic), id(getNextId(lb)), name(std::move(n)), bvSize(bvs),
      cType(type) {}

LogicTerm::LogicTerm(std::string n, const uint64_t identifier, CType type,
                     Logic* logic)
    : lb(logic), id(identifier), name(std::move(n)), cType(type) {}

LogicTerm LogicTerm::noneTerm() {
  return {OpType::None, "None", CType::BOOL, nullptr};
}

LogicTerm LogicTerm::eq(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::EQ, a, b};
}

LogicTerm LogicTerm::neq(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::XOR, a, b};
}

LogicTerm LogicTerm::o(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::OR, a, b};
}

LogicTerm LogicTerm::a(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::AND, a, b};
}

LogicTerm LogicTerm::bvAnd(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::BitAnd, a, b};
}

LogicTerm LogicTerm::bvOr(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::BitOr, a, b};
}

LogicTerm LogicTerm::bvXor(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::BitXor, a, b};
}

LogicTerm LogicTerm::implies(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::IMPL, a, b};
}

LogicTerm LogicTerm::ite(const LogicTerm& a, const LogicTerm& b,
                         const LogicTerm& c) {
  return {OpType::ITE, a, b, c};
}

LogicTerm LogicTerm::add(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::ADD, a, b};
}

LogicTerm LogicTerm::sub(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::SUB, a, b};
}

LogicTerm LogicTerm::mul(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::MUL, a, b};
}

LogicTerm LogicTerm::div(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::DIV, a, b};
}

LogicTerm LogicTerm::gt(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::GT, a, b};
}

LogicTerm LogicTerm::lt(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::LT, a, b};
}

LogicTerm LogicTerm::gte(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::GTE, a, b};
}

LogicTerm LogicTerm::lte(const LogicTerm& a, const LogicTerm& b) {
  return {OpType::LTE, a, b};
}

LogicTerm LogicTerm::neg(const LogicTerm& a) { return {OpType::NEG, a}; }

LogicTerm LogicTerm::operator&&(const LogicTerm& other) const {
  return LogicTerm::a(*this, other);
}

LogicTerm LogicTerm::operator&(const LogicTerm& other) const {
  return LogicTerm::bvAnd(*this, other);
}

LogicTerm LogicTerm::operator|(const LogicTerm& other) const {
  return LogicTerm::bvOr(*this, other);
}

LogicTerm LogicTerm::operator^(const LogicTerm& other) const {
  return LogicTerm::bvXor(*this, other);
}

LogicTerm LogicTerm::operator||(const LogicTerm& other) const {
  return LogicTerm::o(*this, other);
}

LogicTerm LogicTerm::operator==(const LogicTerm& other) const {
  return LogicTerm::eq(*this, other);
}

LogicTerm LogicTerm::operator!=(const LogicTerm& other) const {
  return LogicTerm::neq(*this, other);
}

LogicTerm LogicTerm::operator+(const LogicTerm& other) const {
  return LogicTerm::add(*this, other);
}

LogicTerm LogicTerm::operator-(const LogicTerm& other) const {
  return LogicTerm::sub(*this, other);
}

LogicTerm LogicTerm::operator*(const LogicTerm& other) const {
  return LogicTerm::mul(*this, other);
}

LogicTerm LogicTerm::operator/(const LogicTerm& other) const {
  return LogicTerm::div(*this, other);
}

LogicTerm LogicTerm::operator>(const LogicTerm& other) const {
  return LogicTerm::gt(*this, other);
}

LogicTerm LogicTerm::operator<(const LogicTerm& other) const {
  return LogicTerm::lt(*this, other);
}

LogicTerm LogicTerm::operator>=(const LogicTerm& other) const {
  return LogicTerm::gte(*this, other);
}

LogicTerm LogicTerm::operator<=(const LogicTerm& other) const {
  return LogicTerm::lte(*this, other);
}

LogicTerm LogicTerm::operator!() const { return LogicTerm::neg(*this); }

bool LogicTerm::isConst() const { return getOpType() == OpType::Constant; }

bool LogicTerm::getBoolValue() const {
  switch (cType) {
  case CType::BOOL:
    return value;
  case CType::INT:
    return iValue != 0;
  case CType::REAL:
    return fValue != 0;
  case CType::BITVECTOR:
    return bvValue != 0;
  default:
    return false;
  }
}

int LogicTerm::getIntValue() const {
  switch (cType) {
  case CType::BOOL:
    return value ? 1 : 0;
  case CType::INT:
    return iValue;
  case CType::REAL:
    return static_cast<int32_t>(std::floor(fValue));
  case CType::BITVECTOR:
    return static_cast<int>(bvValue);
  default:
    return std::numeric_limits<int>::max();
  }
}

double LogicTerm::getFloatValue() const {
  switch (cType) {
  case CType::BOOL:
    return value ? 1.0 : 0.0;
  case CType::INT:
    return iValue;
  case CType::REAL:
    return fValue;
  case CType::BITVECTOR:
    return static_cast<double>(bvValue);
  default:
    return std::numeric_limits<double>::max();
  }
}

uint64_t LogicTerm::getBitVectorValue() const {
  switch (cType) {
  case CType::BOOL:
    return value ? 1.0 : 0.0;
  case CType::INT:
    return static_cast<uint64_t>(iValue);
  case CType::REAL:
    return static_cast<uint64_t>(fValue);
  case CType::BITVECTOR:
    return bvValue & (static_cast<uint64_t>(std::pow(2, bvSize)) - 1U);
  default:
    return std::numeric_limits<uint64_t>::max();
  }
}

uint16_t LogicTerm::getBitVectorSize() const {
  switch (cType) {
  case CType::BOOL:
    return 1U;
  case CType::INT:
    return 32U;
  case CType::REAL:
    return 256U;
  case CType::BITVECTOR:
    return bvSize;
  default:
    return std::numeric_limits<uint16_t>::max();
  }
}

bool LogicTerm::deepEquals(const LogicTerm& other) const {
  if (getOpType() == OpType::Variable && getID() == other.getID()) {
    return true;
  }
  if (getDepth() != other.getDepth()) {
    return false;
  }
  if (getOpType() != other.getOpType()) {
    return false;
  }
  if (getName() != other.getName()) {
    return false;
  }
  if (getNodes().size() != other.getNodes().size()) {
    return false;
  }
  if (getCType() != other.getCType()) {
    return false;
  }
  for (size_t i = 0U; i < getNodes().size(); ++i) {
    if (!getNodes()[i].deepEquals(other.getNodes()[i])) {
      return false;
    }
  }
  return true;
}

Logic* LogicTerm::getValidLogicPtr(const LogicTerm& a, const LogicTerm& b) {
  if (a.isConst() || b.isConst()) {
    if (!a.isConst()) {
      return a.getLogic();
    }
    if (!b.isConst()) {
      return b.getLogic();
    }
    return nullptr;
  }
  if (a.getLogic() == b.getLogic()) {
    return a.getLogic();
  }
  throw std::runtime_error("Logic mismatch");
}

Logic* LogicTerm::getValidLogicPtr(const LogicTerm& a, const LogicTerm& b,
                                   const LogicTerm& c) {
  if (a.isConst() || b.isConst() || c.isConst()) {
    if (!a.isConst()) {
      return a.getLogic();
    }
    if (!b.isConst()) {
      return b.getLogic();
    }
    if (!c.isConst()) {
      return c.getLogic();
    }
    return nullptr;
  }
  if (a.getLogic() == b.getLogic() && b.getLogic() == c.getLogic()) {
    return a.getLogic();
  }
  throw std::runtime_error("Logic mismatch");
}

LogicTerm LogicTerm::combineTerms(const LogicTerm& a, const LogicTerm& b,
                                  OpType op, Logic* logic) {
  if ((a.getOpType() == op || b.getOpType() == op) && isAssociative(op)) {
    std::vector<LogicTerm> terms{};
    terms.reserve(a.getNodes().size() + b.getNodes().size());

    auto res = getFlatTerms(a, op);
    terms.insert(terms.end(), res.begin(), res.end());
    res = getFlatTerms(b, op);
    terms.insert(terms.end(), res.begin(), res.end());

    return {op, terms, getTargetCType(a, b, op), logic};
  }
  return {op, a, b, getTargetCType(a, b, op), logic};
}

LogicTerm LogicTerm::combineConst(const LogicTerm& a, const LogicTerm& b,
                                  OpType op, Logic* logic) {
  if (!a.isConst() && !b.isConst()) {
    // erroneous function call
    throw std::runtime_error("Both terms are not constants");
  }
  if (a.isConst() && b.isConst()) {
    // combine the values, return new const
    switch (op) {
    case OpType::AND:
      return LogicTerm(a.getBoolValue() && b.getBoolValue());
    case OpType::OR:
      return LogicTerm(a.getBoolValue() || b.getBoolValue());
    case OpType::IMPL:
      return LogicTerm(!a.getBoolValue() || b.getBoolValue());
    case OpType::ADD:
      return LogicTerm(a.getIntValue() + b.getIntValue());
    case OpType::SUB:
      return LogicTerm(a.getIntValue() - b.getIntValue());
    case OpType::MUL:
      return LogicTerm(a.getIntValue() * b.getIntValue());
    case OpType::DIV:
      // NOLINTNEXTLINE(clang-analyzer-core.DivideZero)
      return LogicTerm(a.getIntValue() / b.getIntValue());
    case OpType::GT:
      return LogicTerm(a.getIntValue() > b.getIntValue());
    case OpType::LT:
      return LogicTerm(a.getIntValue() < b.getIntValue());
    case OpType::GTE:
      return LogicTerm(a.getIntValue() >= b.getIntValue());
    case OpType::LTE:
      return LogicTerm(a.getIntValue() <= b.getIntValue());
    case OpType::EQ:
      return LogicTerm(a.getFloatValue() == b.getFloatValue());
    case OpType::XOR:
      return LogicTerm(a.getFloatValue() != b.getFloatValue());
    default:
      throw std::runtime_error("Invalid operator");
    }
  } else if (a.isConst() && isCommutative(op)) {
    // since the combineOneConst ignores order of operands
    return combineOneConst(a, b, op, logic);
  } else if (b.isConst()) {
    // const comes at the end anyway
    return combineOneConst(b, a, op, logic);
  } else {
    // combineTerms respects order of operands
    return combineTerms(a, b, op, logic);
  }
}

LogicTerm LogicTerm::combineOneConst(const LogicTerm& constant,
                                     const LogicTerm& other, OpType op,
                                     Logic* logic) {
  switch (op) { // TODO handle other CTypes
  case OpType::AND: {
    if (constant.getBoolValue()) {
      return other;
    }
    return LogicTerm(false);
  }
  case OpType::OR: {
    if (!constant.getBoolValue()) {
      return other;
    }
    return LogicTerm(true);
  }
  case OpType::ADD: {
    if (constant.getFloatValue() == 0) {
      return other;
    }
    return {OpType::ADD, other, constant, CType::INT, logic};
  }
  case OpType::SUB: {
    if (constant.getFloatValue() == 0) {
      return other;
    }
    return {OpType::SUB, other, constant, CType::INT, logic};
  }
  case OpType::MUL: {
    if (constant.getFloatValue() == 0) {
      return LogicTerm(0);
    }
    if (constant.getFloatValue() == 1) {
      return other;
    }
    return {OpType::MUL, other, constant, CType::INT, logic};
  }
  case OpType::DIV: {
    if (constant.getFloatValue() == 0) {
      throw std::runtime_error("Divide by zero");
    }
    if (constant.getFloatValue() == 1) {
      return other;
    }
    return {OpType::DIV, other, constant, CType::INT, logic};
  }
  default: // TODO there are multiple ctypes
    return {op, other, constant, getResultCType(op), logic};
  }
}

LogicTerm LogicTerm::negateTerm(const LogicTerm& term, Logic* logic) {
  if (term.isConst()) {
    switch (term.getCType()) {
    case CType::BOOL:
      return LogicTerm(!term.getBoolValue());
    case CType::INT:
      return LogicTerm(-term.getIntValue());
    case CType::REAL:
      return LogicTerm(-term.getFloatValue());
    default:
      throw std::runtime_error("Invalid CType");
    }
  } else if (term.getOpType() == OpType::NEG) {
    return term.getNodes().front();
  } else {
    return {OpType::NEG, term, term.getCType(), logic};
  };
}

LogicTerm LogicTerm::getBoolConversion() const {
  if (isConst()) {
    if (getCType() == CType::BOOL) {
      return *this;
    }
    return LogicTerm(getBoolValue());
  }
  {
    switch (getCType()) {
    case CType::BOOL:
      return *this;
    case CType::INT:
    case CType::REAL:
      return eq(*this, LogicTerm(0));
    case CType::BITVECTOR:
      return eq(*this, LogicTerm(0, getBitVectorSize()));
    default:
      throw std::runtime_error("Invalid CType");
    }
  }
}

CType LogicTerm::getTargetCType(const LogicTerm& a, const LogicTerm& b,
                                OpType op) {
  if (op == OpType::EQ || op == OpType::XOR || op == OpType::AND ||
      op == OpType::OR || op == OpType::GT || op == OpType::LT ||
      op == OpType::GTE || op == OpType::LTE) {
    return CType::BOOL;
  }
  if (a.getCType() == CType::REAL || b.getCType() == CType::REAL) {
    return CType::REAL;
  }
  if (a.getCType() == CType::BITVECTOR || b.getCType() == CType::BITVECTOR) {
    return CType::BITVECTOR;
  }
  if (a.getCType() == CType::INT || b.getCType() == CType::INT) {
    return CType::INT;
  }
  return CType::BOOL;
}

CType LogicTerm::getTargetCType(CType targetType, const LogicTerm& b) {
  if (targetType == CType::REAL || b.getCType() == CType::REAL) {
    return CType::REAL;
  }
  if (targetType == CType::BITVECTOR || b.getCType() == CType::BITVECTOR) {
    return CType::BITVECTOR;
  }
  if (targetType == CType::INT || b.getCType() == CType::INT) {
    return CType::INT;
  }
  return CType::BOOL;
}

std::vector<LogicTerm> LogicTerm::getFlatTerms(const LogicTerm& t, OpType op) {
  std::vector<LogicTerm> terms;
  if (t.getOpType() != op) {
    terms.push_back(t);
  } else {
    for (const LogicTerm& it : t.getNodes()) {
      if (it.getOpType() != op) {
        terms.push_back(it);
      } else {
        auto res = getFlatTerms(it, op);
        terms.insert(terms.end(), res.begin(), res.end());
      }
    }
  }
  return terms;
}

uint64_t LogicTerm::getMax(const std::vector<LogicTerm>& terms) {
  uint64_t ret = 0;
  for (const auto& it : terms) {
    ret = std::max(ret, it.getDepth());
  }
  return ret + 1;
}

uint16_t LogicTerm::getMaxBVSize(const std::vector<LogicTerm>& terms) {
  uint16_t ret = 0U;
  for (const auto& it : terms) {
    ret = std::max(ret, it.getBitVectorSize());
  }
  return ret;
}

std::size_t TermHash::operator()(const LogicTerm& t) const {
  if (t.getOpType() == OpType::None) {
    throw std::runtime_error("Invalid OpType");
  }
  return t.getID();
}

bool TermHash::operator()(const LogicTerm& t1, const LogicTerm& t2) const {
  if (t1.getOpType() == OpType::None || t2.getOpType() == OpType::None) {
    throw std::runtime_error("Invalid OpType");
  }
  if (t1.getOpType() == OpType::Constant &&
      t2.getOpType() == OpType::Constant && t1.getCType() == t2.getCType()) {
    switch (t1.getCType()) {
    case CType::BOOL:
      return t1.getBoolValue() == t2.getBoolValue();
    case CType::INT:
      return t1.getIntValue() == t2.getIntValue();
    case CType::REAL:
      return t1.getFloatValue() == t2.getFloatValue();
    case CType::BITVECTOR:
      return t1.getBitVectorValue() == t2.getBitVectorValue();
    default:
      return false;
    }
  } else {
    return t1.getID() == t2.getID();
  }
}

bool TermDepthComparator::operator()(const LogicTerm& t1,
                                     const LogicTerm& t2) const {
  if (t1.getOpType() == OpType::None || t2.getOpType() == OpType::None) {
    throw std::runtime_error("Invalid OpType");
  }
  if ((t1.getOpType() == OpType::Variable &&
       t2.getOpType() == OpType::Variable) ||
      (t1.getOpType() == OpType::Constant &&
       t2.getOpType() == OpType::Constant) ||
      t1.getDepth() == t2.getDepth()) {
    return t1.getID() > t2.getID();
  }
  {
    return t1.getDepth() > t2.getDepth();
  }
}
} // namespace logicbase
