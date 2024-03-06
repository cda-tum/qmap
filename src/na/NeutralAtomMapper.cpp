//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "NeutralAtomMapper.hpp"

#include "Architecture.hpp"
#include "Graph.hpp"
#include "Layer.hpp"
#include "QuantumComputation.hpp"
#include "na/Definitions.hpp"
#include "operations/NAGlobalOperation.hpp"
#include "operations/NALocalOperation.hpp"
#include "operations/NAShuttlingOperation.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include <__algorithm/remove.h>
#include <__algorithm/remove_if.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

namespace na {

auto NeutralAtomMapper::preprocess() -> void {
  // validate circuit
  for (const auto& op : initialQc) {
    if (op->isStandardOperation()) {
      if (!isIndividual(*op) && op->getNcontrols() + op->getNtargets() > 2) {
        throw std::logic_error("Operations acting on more than two qubits "
                               "are not supported yet.");
      }
      if (!arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
        if (!arch.isAllowedGlobally({op->getType(), op->getNcontrols()})) {
          std::stringstream ss;
          ss << "The chosen architecture does not support the operation "
             << OpType{op->getType(), 0} << " either locally or globally.";
          throw std::invalid_argument(ss.str());
        }
        // the gate is global: if it is a 1Q-gate it must be applied globally
        if (isIndividual(*op) and !isGlobal(*op)) {
          std::stringstream ss;
          ss << "The chosen architecture does not support the operation "
             << OpType{op->getType(), op->getNcontrols()} << " locally.";
          throw std::invalid_argument(ss.str());
        }
      }
    } else {
      throw std::logic_error("Operation type is not supported.");
    }
  }
}

auto NeutralAtomMapper::postprocess() -> void {
  // make patches
  const auto logicQC = mappedQc;
  mappedQc.clear();
  mappedQc.clearInitialPositions();
  const auto rows = static_cast<std::int64_t>(config.getPatchRows());
  const auto cols = static_cast<std::int64_t>(config.getPatchCols());
  for (const auto& p : logicQC.getInitialPositions()) {
    for (std::int64_t r = 0; r < rows; ++r) {
      for (std::int64_t c = 0; c < cols; ++c) {
        mappedQc.emplaceInitialPosition(
            std::make_shared<Point>(initialArch.getPositionOffsetBy(*p, r, c)));
      }
    }
  }
  for (const auto& op : logicQC) {
    if (op->isGlobalOperation()) {
      mappedQc.emplaceBack(op->clone());
    } else if (op->isLocalOperation()) {
      const auto& lop = dynamic_cast<NALocalOperation&>(*op);
      for (std::int64_t r = 0; r < rows; ++r) {
        std::vector<std::shared_ptr<Point>> positions;
        for (const auto& p : lop.getPositions()) {
          for (std::int64_t c = 0; c < cols; ++c) {
            positions.emplace_back(std::make_shared<Point>(
                initialArch.getPositionOffsetBy(*p, r, c)));
          }
        }
        mappedQc.emplaceBack<NALocalOperation>(lop.getType(), lop.getParams(),
                                               positions);
      }
    } else if (op->isShuttlingOperation()) {
      const auto& sop = dynamic_cast<NAShuttlingOperation&>(*op);
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      for (std::size_t i = 0; i < sop.getStart().size(); ++i) {
        const auto s = sop.getStart()[i];
        const auto e = sop.getEnd()[i];
        for (std::int64_t r = 0; r < rows; ++r) {
          for (std::int64_t c = 0; c < cols; ++c) {
            start.emplace_back(std::make_shared<Point>(
                initialArch.getPositionOffsetBy(*s, r, c)));
            end.emplace_back(std::make_shared<Point>(
                initialArch.getPositionOffsetBy(*e, r, c)));
          }
        }
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(sop.getType(), start, end);
    } else {
      throw std::logic_error("Operation type is not supported.");
    }
  }
  const auto prelQC = mappedQc;
  mappedQc.clear();
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  for (const auto& op : prelQC) {
    if (op->isShuttlingOperation()) {
      const auto& sop = dynamic_cast<NAShuttlingOperation&>(*op);
      if (sop.getType() == MOVE) {
        std::vector<std::shared_ptr<Point>> hOffsetStart;
        std::vector<std::shared_ptr<Point>> hOffsetEnd;
        std::vector<std::shared_ptr<Point>> vMoveStart;
        std::vector<std::shared_ptr<Point>> vMoveEnd;
        std::vector<std::shared_ptr<Point>> hMoveStart;
        std::vector<std::shared_ptr<Point>> hMoveEnd;
        std::vector<std::shared_ptr<Point>> vOffsetStart;
        std::vector<std::shared_ptr<Point>> vOffsetEnd;
        bool                                vOffset = false;
        for (std::size_t i = 0; i < sop.getStart().size(); ++i) {
          auto        start = *sop.getStart()[i];
          auto        end   = *sop.getEnd()[i];
          const auto  dx    = end.x - start.x;
          Point const mid   = {start.x, end.y};
          if (dx > 0) {
            try {
              const auto& s = arch.getNearestSiteRight(mid, true);
              if (arch.getPositionOfSite(s).x < end.x) {
                vOffset = true;
              }
            } catch (std::invalid_argument& e) {
            }
          } else if (dx < 0) {
            try {
              const auto& s = arch.getNearestSiteLeft(mid, true);
              if (arch.getPositionOfSite(s).x > end.x) {
                vOffset = true;
              }
            } catch (std::invalid_argument& e) {
            }
          }
        }
        for (std::size_t i = 0; i < sop.getStart().size(); ++i) {
          auto       start = *sop.getStart()[i];
          auto       end   = *sop.getEnd()[i];
          const auto dx    = end.x - start.x;
          const auto dy    = end.y - start.y;
          if (dy > 0) {
            try {
              const auto& s = arch.getNearestSiteDown(start, true);
              if (arch.getPositionOfSite(s).y < end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Point>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Point>(start));
              }
            } catch (std::invalid_argument& e) {
            }
          } else if (dy < 0) {
            try {
              const auto& s = arch.getNearestSiteUp(start, true);
              if (arch.getPositionOfSite(s).y > end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Point>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Point>(start));
              }
            } catch (std::invalid_argument& e) {
            }
          }
          Point mid = {start.x, end.y};
          if (vOffset) {
            mid.y += (dy >= 0 ? -d : d);
          }
          if (start.y != mid.y) {
            vMoveStart.emplace_back(std::make_shared<Point>(start));
            start = mid;
            vMoveEnd.emplace_back(std::make_shared<Point>(start));
          }
          if (start.x != end.x) {
            hMoveStart.emplace_back(std::make_shared<Point>(start));
            start.x = end.x;
            hMoveEnd.emplace_back(std::make_shared<Point>(start));
          }
          if (start.y != end.y) {
            vOffsetStart.emplace_back(std::make_shared<Point>(start));
            start.y = end.y;
            vOffsetEnd.emplace_back(std::make_shared<Point>(start));
          }
        }
        if (!hOffsetStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, hOffsetStart,
                                                     hOffsetEnd);
        }
        if (!vMoveStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, vMoveStart,
                                                     vMoveEnd);
        }
        if (!hMoveStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, hMoveStart,
                                                     hMoveEnd);
        }
        if (!vOffsetStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, vOffsetStart,
                                                     vOffsetEnd);
        }
      } else {
        mappedQc.emplaceBack(op->clone());
      }
    } else {
      mappedQc.emplaceBack(op->clone());
    }
  }
}

auto NeutralAtomMapper::checkApplicability(
    const std::unique_ptr<qc::Operation>& op,
    const std::vector<Atom>&              placement) const -> bool {
  assert(op->isStandardOperation()); // ensured by preprocess
  if (arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
    if (op->getNcontrols() == 0) {
      // individual gate that can act on one or more atoms
      return std::all_of(
          op->getTargets().cbegin(), op->getTargets().cend(),
          [&](const auto& qubit) {
            switch (placement[qubit].positionStatus) {
            case Atom::PositionStatus::UNDEFINED:
              // check whether the gate is applicable in one of the currently
              // selected zones
              return std::any_of(
                  placement[qubit].zones.cbegin(),
                  placement[qubit].zones.cend(), [&](const auto& z) {
                    return arch.isAllowedLocally({op->getType(), 0}, z);
                  });
            case Atom::PositionStatus::DEFINED:
              // check whether the gate is applicable at the current position
              return arch.isAllowedLocallyAt({op->getType(), 0},
                                             *placement[qubit].currentPosition);
            }
          });
    } else {
      // TODO single controlled local gate that acts exactly on two atoms
      return false;
    }
  } else {
    assert(arch.isAllowedGlobally({op->getType(), op->getNcontrols()}));
    if (isIndividual(*op)) {
      assert(isGlobal(*op));
      // purely globally gate can always be applied
      return true;
    } else {
      assert(op->getNcontrols() + op->getNtargets() == 2);
      // TODO single controlled global gate that acts exactly on two atoms
      return false;
    }
  }
}

auto NeutralAtomMapper::updatePlacement(
    const std::unique_ptr<qc::Operation>& op,
    std::vector<Atom>&                    placement) const -> void {
  if (arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
    if (op->getNcontrols() == 0) {
      // individual gate that can act on one or more atoms
      std::for_each(
          op->getTargets().cbegin(), op->getTargets().cend(),
          [&](const auto& qubit) {
            switch (placement[qubit].positionStatus) {
            case Atom::PositionStatus::UNDEFINED:
              // remove all zones where the gate is not applicable
              placement[qubit].zones.erase(
                  std::remove_if(
                      placement[qubit].zones.begin(),
                      placement[qubit].zones.end(),
                      [&](auto& z) {
                        return !arch.isAllowedLocally({op->getType(), 0}, z);
                      }),
                  placement[qubit].zones.end());
              break;
            case Atom::PositionStatus::DEFINED:
              break;
            }
          });
    } else {
      // TODO single controlled local gate that acts exactly on two atoms
    }
  } else {
    assert(arch.isAllowedGlobally({op->getType(), op->getNcontrols()}));
    if (isIndividual(*op)) {
      assert(isGlobal(*op));
      // purely globally gate can always be applied
      std::for_each(
          op->getTargets().cbegin(), op->getTargets().cend(),
          [&](const auto qubit) {
            switch (placement[qubit].positionStatus) {
            case Atom::PositionStatus::UNDEFINED:
              placement[qubit].zones.erase(
                  std::remove_if(
                      placement[qubit].zones.begin(),
                      placement[qubit].zones.end(),
                      [&](auto& z) {
                        return !arch.isAllowedGlobally({op->getType(), 0}, z);
                      }),
                  placement[qubit].zones.end());
              break;
            case Atom::PositionStatus::DEFINED:
              break;
            }
          });
    } else {
      assert(op->getNcontrols() + op->getNtargets() == 2);
      // TODO single controlled global gate that acts exactly on two atoms
    }
  }
}

[[nodiscard]] auto NeutralAtomMapper::getMisplacement(
    const std::vector<Atom>&                           initial,
    const std::unordered_map<qc::Qubit, std::int64_t>& target,
    const qc::Qubit& q) const -> std::int64_t {
  if (initial.at(q).positionStatus == Atom::PositionStatus::DEFINED) {
    return std::accumulate(target.cbegin(), target.cend(), 0,
                           [&](const auto& acc, const auto& p) {
                             if (initial.at(p.first).positionStatus ==
                                 Atom::PositionStatus::DEFINED) {
                               if (initial.at(p.first).currentPosition->x >
                                       initial.at(q).currentPosition->x and
                                   p.second < target.at(q)) {
                                 return acc + 1;
                               }
                               if (initial.at(p.first).currentPosition->x <
                                       initial.at(q).currentPosition->x and
                                   p.second > target.at(q)) {
                                 return acc - 1;
                               }
                             }
                             return acc;
                           }) +
           // count the spots the atom must be moved relatively
           std::count_if(
               target.cbegin(), target.cend(),
               [&](const auto& p) { return p.second < target.at(q); }) -
           std::count_if(target.cbegin(), target.cend(), [&](const auto& p) {
             return initial[p.first].currentPosition->x <
                    initial[q].currentPosition->x;
           });
  }
  return 0;
}

auto NeutralAtomMapper::map(const qc::QuantumComputation& qc) -> void {
  auto startPreprocess = std::chrono::high_resolution_clock::now();
  initialQc    = qc;
  auto nqubits = initialQc.getNqubits();
  mappedQc     = NAQuantumComputation();
  preprocess();
  // store the placement of atoms, both the initial one (needed later) and the
  // current one leave atoms unmapped as long as possible. This mighty induce
  // some restrictions on the mapping in which zone the atom may be placed.
  const std::vector<Zone> initialZones = arch.getInitialZones();
  std::vector<Atom>       placement(nqubits);
  std::for_each(placement.begin(), placement.end(),
                [&](auto& p) { p = Atom(initialZones); });
  std::vector<bool>             initialFreeSites(arch.getNSites(), true);
  std::vector<bool>             currentFreeSites(arch.getNSites(), true);
  std::unordered_set<qc::Qubit> currentlyShuttling{};
  // get start time
  auto startMapping = std::chrono::high_resolution_clock::now();
  //============================ START MAPPING ============================
  Layer const layer         = Layer(initialQc);
  const auto& executableSet = layer.getExecutableSet();
  auto        it            = (*executableSet)->begin();
  while (it != (*executableSet)->end()) {
    // 1. execute all gates that are directly applicable and do not need
    //    shuttling
    while (it != (*executableSet)->end()) {
      const std::unique_ptr<qc::Operation>& op = *(*it)->getOperation();
      if (checkApplicability(op, placement)) {
        updatePlacement(op, placement);
        (*it)->execute();
        if (isGlobal(*op) and
            arch.isAllowedGlobally({op->getType(), op->getNcontrols()})) {
          mappedQc.emplaceBack<NAGlobalOperation>(
              OpType{op->getType(), op->getNcontrols()}, op->getParameter());
        } else {
          // collect executable gates of the same type
          std::vector<std::shared_ptr<Point>> positions = {
              placement[op->getTargets().front()].currentPosition};
          for (const auto& v : layer.getExecutablesOfType(
                   {op->getType(), op->getNcontrols()})) {
            const auto& op2 = *v->getOperation();
            if (checkApplicability(op2, placement) and
                op->getParameter() == op2->getParameter()) {
              updatePlacement(op2, placement);
              v->execute();
              positions.emplace_back(
                  placement[op2->getTargets().front()].currentPosition);
            }
          }
          mappedQc.emplaceBack<NALocalOperation>(
              OpType{op->getType(), op->getNcontrols()}, op->getParameter(),
              positions);
        }
        it = (*executableSet)->begin();
      } else {
        it++;
      }
    }
    it = (*executableSet)->begin();
    if (it == (*executableSet)->end()) {
      break;
    }
    // 2. when no such gates are left, extract an interaction graph of gates
    //    of the same type and two targets, i.e. cz gates
    const Graph<std::shared_ptr<Layer::DAGVertex>> graph =
        layer.constructInteractionGraph({qc::OpType::Z, 1});
    if (graph.getNVertices() == 0) {
      throw std::logic_error(
          "Other gates than cz are not supported for mapping yet.");
      // TODO: support other gates than cz
    }
    const auto& sequence = graph.computeSequence();
    const auto& moveable = sequence.first;
    const auto& fixed    = sequence.second;
    // 3. move the atoms accordingly and execute the gates
    // this distance is used for spacing atoms that should interact or pass
    // another atom
    const auto d  = static_cast<std::int64_t>(arch.getMinAtomDistance());
    const auto dx = static_cast<std::int64_t>(config.getPatchCols()) * static_cast<std::int64_t>(arch.getNoInteractionRadius());
    // pick up the fixed atoms and move them to the interaction zone
    // get a vector of the fixed atoms ordered by their initial position from
    // left to right
    std::vector<qc::Qubit> fixedOrdered;
    std::transform(fixed.cbegin(), fixed.cend(), back_inserter(fixedOrdered),
                   [](const auto& v) { return v.first; });
    std::sort(fixedOrdered.begin(), fixedOrdered.end(),
              [&](const auto& a, const auto& b) {
                return fixed.at(a) < fixed.at(b);
              });
    // get a vector of the fixed atoms in the order to pick them up based on
    // their misplacement value
    std::vector<qc::Qubit> pickUpOrderFixed(fixedOrdered);
    std::sort(pickUpOrderFixed.begin(), pickUpOrderFixed.end(),
              [&](const auto& a, const auto& b) {
                return std::abs(getMisplacement(placement, fixed, a)) >
                       std::abs(getMisplacement(placement, fixed, b));
              });
    // repeat until all atoms are picked up
    while (!pickUpOrderFixed.empty()) {
      // extract the first atom from the pick up order
      const auto q = pickUpOrderFixed.front();
      pickUpOrderFixed.erase(
          std::remove(pickUpOrderFixed.begin(), pickUpOrderFixed.end(), q),
          pickUpOrderFixed.end());
      // get the index of the picked atom among the fixed atoms ordered by their
      // initial position
      const auto& fixIt =
          std::find(fixedOrdered.cbegin(), fixedOrdered.cend(), q);
      const auto& i = static_cast<std::uint64_t>(
          std::distance(fixedOrdered.cbegin(), fixIt));
      // if the placement of the atom is undefined, find a good placement for it
      if (placement[q].positionStatus == Atom::PositionStatus::UNDEFINED) {
        // calculate not picked up atoms to the left in the resulting order
        // note: all remaining atoms that are not picked up yet have undefined
        // positions as well
        auto notPickedUpLeft = 0UL;
        for (std::size_t j = 0; j < i; ++j) {
          const auto& p = fixedOrdered[j];
          // not picked up yet and left
          // of q in the end
          if (currentlyShuttling.find(p) == currentlyShuttling.cend() and
              fixed.at(p) < fixed.at(q)) {
            ++notPickedUpLeft;
          }
        }
        Zone        zone           = 0;
        Index       row            = 0;
        std::size_t freeSpotsInRow = 0;
        for (const auto& z : placement[q].zones) {
          for (std::size_t r = 0; r < arch.getNrowsInZone(z); ++r) {
            const auto& sitesInRow = arch.getSitesInRow(z, r);
            const auto& n =
                std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                                [&](const auto& acc, const auto& s) {
                                  return acc + (initialFreeSites[s] ? 1 : 0);
                                });
            if (freeSpotsInRow <= notPickedUpLeft and n > freeSpotsInRow) {
              freeSpotsInRow = n;
              zone           = z;
              row            = r;
            }
          }
        }
        // place q on the i-th free spot where i is the minimum of the number of
        // not picked up atoms to the left and free spots available
        auto possibleSites = arch.getSitesInRow(zone, row);
        possibleSites.erase(
            std::remove_if(possibleSites.begin(), possibleSites.end(),
                           [&](const auto s) { return !initialFreeSites[s]; }),
            possibleSites.end());
        const auto s =
            possibleSites[std::min(notPickedUpLeft, freeSpotsInRow - 1)];
        placement[q].positionStatus   = Atom::PositionStatus::DEFINED;
        *placement[q].initialPosition = arch.getPositionOfSite(s);
        initialFreeSites[s]           = false;
        currentFreeSites[s]           = false;
      }
      // here the position of q is defined
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      std::vector<std::shared_ptr<Point>> loadStart;
      std::vector<std::shared_ptr<Point>> loadEnd;
      const auto currentX = placement[q].currentPosition->x;
      const auto y        = placement[q].currentPosition->y;
      // pick up q itself
      loadStart.emplace_back(placement[q].currentPosition);
      currentFreeSites[arch.getSiteAt(*placement[q].currentPosition)] = true;
      placement[q].currentPosition = std::make_shared<Point>(currentX + d, y);
      loadEnd.emplace_back(placement[q].currentPosition);
      currentlyShuttling.insert(q);
      // iterate through all not yet picked up atoms to the left and check
      // whether they can be picked up
      auto x = currentX - d;
      if (i > 0) {
        for (std::size_t j = i - 1;; --j) {
          const auto p = fixedOrdered[j];
          // if the atom is picked up,
          // move it to the correct row
          if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
            start.emplace_back(placement[p].currentPosition);
            placement[p].currentPosition = std::make_shared<Point>(x, y);
            end.emplace_back(placement[p].currentPosition);
            const auto nx = arch.getNearestXLeft(x);
            x             = nx == x ? x - dx : nx - d;
          } else {
            // check whether j can be
            // picked up
            if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
              if (placement[p].currentPosition->y == y and
                  placement[p].currentPosition->x <= x - d) {
                // pick up p
                pickUpOrderFixed.erase(std::remove(pickUpOrderFixed.begin(),
                                                   pickUpOrderFixed.end(), p),
                                       pickUpOrderFixed.end());
                x = placement[p].currentPosition->x;
                loadStart.emplace_back(placement[p].currentPosition);
                currentFreeSites[arch.getSiteAt(
                    *placement[p].currentPosition)] = true;
                placement[p].currentPosition =
                    std::make_shared<Point>(x + d, y);
                loadEnd.emplace_back(placement[p].currentPosition);
                currentlyShuttling.insert(p);
                x -= d;
              }
            } else {
              // if atom is not placed
              // yet find a site to the
              // left that it can be
              // picked up together with
              // q, i.e. find next free
              // site to the left
              std::int64_t freeX = x;
              bool         free  = false;
              try {
                while (!free) {
                  const auto& site = arch.getNearestSiteLeft({freeX, y});
                  freeX            = arch.getPositionOfSite(site).x;
                  if (initialFreeSites[site] and
                      std::find(placement[p].zones.cbegin(),
                                placement[p].zones.cend(),
                                arch.getZoneOfSite(site)) !=
                          placement[p].zones.cend()) {
                    // the site is free
                    // and satisfies the
                    // zone restrictions
                    // of the atom
                    free = true;
                  } else {
                    freeX -= d;
                  }
                }
              } catch (std::invalid_argument& e) {
                // if x reached the left
                // end and there was no
                // site to the left
                // anymore, free remains
                // false in this case
              }
              if (free) {
                // place p on the free
                // site
                placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
                *placement[p].initialPosition = {freeX, y};
                initialFreeSites[arch.getSiteAt(
                    *placement[p].initialPosition)] = false;
                // pick up p
                pickUpOrderFixed.erase(std::remove(pickUpOrderFixed.begin(),
                                                   pickUpOrderFixed.end(), p),
                                       pickUpOrderFixed.end());
                loadStart.emplace_back(placement[p].currentPosition);
                placement[p].currentPosition =
                    std::make_shared<Point>(freeX + d, y);
                loadEnd.emplace_back(placement[p].currentPosition);
                currentlyShuttling.insert(p);
                x = freeX - d;
              }
            }
          }
          if (j == 0) {
            break;
          }
        }
      }
      // iterate through all not yet picked up atoms to the right and check
      // whether they can be picked up
      x = arch.getNearestXRight(currentX + d) + d;
      for (std::size_t j = i + 1; j < fixedOrdered.size(); ++j) {
        const auto p = fixedOrdered[j];
        // if the atom is picked up, move it to the correct row
        if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
          start.emplace_back(placement[p].currentPosition);
          placement[p].currentPosition = std::make_shared<Point>(x, y);
          end.emplace_back(placement[p].currentPosition);
          const auto nx = arch.getNearestXRight(x) + d;
          x             = nx == x ? x + dx : nx + d;
        } else {
          // check whether j can be
          // picked up
          if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
            if (placement[p].currentPosition->y == y and
                placement[p].currentPosition->x >= x - d) {
              // pick up p
              pickUpOrderFixed.erase(std::remove(pickUpOrderFixed.begin(),
                                                 pickUpOrderFixed.end(), p),
                                     pickUpOrderFixed.end());
              x = placement[p].currentPosition->x;
              loadStart.emplace_back(placement[p].currentPosition);
              currentFreeSites[arch.getSiteAt(*placement[p].currentPosition)] =
                  true;
              placement[p].currentPosition = std::make_shared<Point>(x + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
              currentlyShuttling.insert(p);
              const auto nx = arch.getNearestXRight(x + d);
              x             = nx == x + d ? nx - d + dx : nx + d;
            }
          } else {
            // if atom is not placed yet
            // find a site to the right
            // that it can be picked up
            // together with q, i.e. find
            // next free site to the
            // right
            std::int64_t freeX = x;
            bool         free  = false;
            try {
              while (!free) {
                const auto& site = arch.getNearestSiteRight({freeX - d, y});
                freeX            = arch.getPositionOfSite(site).x;
                if (initialFreeSites[site] and
                    std::find(placement[p].zones.cbegin(),
                              placement[p].zones.cend(),
                              arch.getZoneOfSite(site)) !=
                        placement[p].zones.cend()) {
                  // the site is free and
                  // satisfies the zone
                  // restrictions of the
                  // atom
                  free = true;
                } else {
                  freeX = arch.getNearestXRight(freeX + d) + d;
                }
              }
            } catch (std::invalid_argument& e) {
              // if x reached the right
              // end and there was no
              // site to the right
              // anymore, free remains
              // false in this case
            }
            if (free) {
              // place p on the free site
              placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
              *placement[p].initialPosition = {freeX, y};
              initialFreeSites[arch.getSiteAt(*placement[p].initialPosition)] =
                  false;
              // pick up p
              pickUpOrderFixed.erase(std::remove(pickUpOrderFixed.begin(),
                                                 pickUpOrderFixed.end(), p),
                                     pickUpOrderFixed.end());
              loadStart.emplace_back(placement[p].currentPosition);
              placement[p].currentPosition =
                  std::make_shared<Point>(freeX + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
              currentlyShuttling.insert(p);
              const auto nx = arch.getNearestXRight(freeX + d);
              x             = nx == freeX + d ? nx - d + dx : nx + d;
            }
          }
        }
      }
      if (!start.empty()) {
        mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(LOAD, loadStart, loadEnd);
    }
    // all atoms are picked up in order, move them to the interaction zone and
    // store them there
    std::vector<std::shared_ptr<Point>> startFixed;
    std::vector<std::shared_ptr<Point>> midFixed;
    std::vector<std::shared_ptr<Point>> endFixed;
    const Zone                          interactionZone =
        *arch.getPropertiesOfOperation({qc::OpType::Z, 1}).zones.begin();
    const auto sites = arch.getSitesInRow(interactionZone, 0);
    for (const auto& [q, x] : fixed) {
      if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
        std::stringstream ss;
        ss << "Atom " << q << " was unexpectedly not picked up.";
        throw std::logic_error(ss.str());
      }
      if (!currentFreeSites[sites[static_cast<std::size_t>(x)]]) {
        throw std::logic_error(
            "Target site in interaction zone is unexpectedly occupied.");
      }
      const auto pos =
          arch.getPositionOfSite(sites[static_cast<std::size_t>(x)]);
      startFixed.emplace_back(placement[q].currentPosition);
      placement[q].currentPosition = std::make_shared<Point>(pos.x + d, pos.y);
      midFixed.emplace_back(placement[q].currentPosition);
      placement[q].currentPosition = std::make_shared<Point>(pos);
      endFixed.emplace_back(placement[q].currentPosition);
      currentFreeSites[sites[static_cast<std::size_t>(x)]] = false;
    }
    currentlyShuttling.clear();
    mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, startFixed, midFixed);
    mappedQc.emplaceBack<NAShuttlingOperation>(STORE, midFixed, endFixed);
    // -----------------------------------------------------------------
    // pick-up the moveable atoms and move them to the interaction zone
    std::vector<qc::Qubit> moveableOrdered;
    std::transform(moveable[0].cbegin(), moveable[0].cend(),
                   back_inserter(moveableOrdered),
                   [](const auto& v) { return v.first; });
    std::sort(moveableOrdered.begin(), moveableOrdered.end(),
              [&](const auto& a, const auto& b) {
                return moveable[0].at(a) < moveable[0].at(b);
              });
    // get a vector of the fixed atoms in the order to pick them up based on
    // their misplacement value
    std::vector<qc::Qubit> pickUpOrderMoveable(moveableOrdered);
    std::sort(pickUpOrderMoveable.begin(), pickUpOrderMoveable.end(),
              [&](const auto& a, const auto& b) {
                return std::abs(getMisplacement(placement, moveable[0], a)) >
                       std::abs(getMisplacement(placement, moveable[0], b));
              });
    // repeat until all atoms are picked up
    while (!pickUpOrderMoveable.empty()) {
      // extract the first atom from the pick up order
      const auto q = pickUpOrderMoveable.front();
      pickUpOrderMoveable.erase(std::remove(pickUpOrderMoveable.begin(),
                                            pickUpOrderMoveable.end(), q),
                                pickUpOrderMoveable.end());
      // get the index of the picked atom among the fixed atoms ordered by their
      // initial position
      const auto& fixIt =
          std::find(moveableOrdered.cbegin(), moveableOrdered.cend(), q);
      const auto& i = static_cast<std::uint64_t>(
          std::distance(moveableOrdered.cbegin(), fixIt));
      // if the placement of the atom is undefined, find a good placement for it
      if (placement[q].positionStatus == Atom::PositionStatus::UNDEFINED) {
        // calculate not picked up atoms to the left in the resulting order
        // note: all remaining atoms that are not picked up yet have undefined
        // positions as well
        auto notPickedUpLeft = 0UL;
        for (std::size_t j = 0; j < i; ++j) {
          const auto& p = moveableOrdered[j];
          // not picked up yet and left
          // of q in the end
          if (currentlyShuttling.find(p) == currentlyShuttling.cend() and
              fixed.at(p) < fixed.at(q)) {
            ++notPickedUpLeft;
          }
        }
        Zone        zone           = 0;
        Index       row            = 0;
        std::size_t freeSpotsInRow = 0;
        for (const auto& z : placement[q].zones) {
          for (std::size_t r = 0; r < arch.getNrowsInZone(z); ++r) {
            const auto& sitesInRow = arch.getSitesInRow(z, r);
            const auto& n =
                std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                                [&](const auto& acc, const auto& s) {
                                  return acc + (initialFreeSites[s] ? 1 : 0);
                                });
            if (freeSpotsInRow <= notPickedUpLeft and n > freeSpotsInRow) {
              freeSpotsInRow = n;
              zone           = z;
              row            = r;
            }
          }
        }
        // place q on the i-th free spot where i is the minimum of the number of
        // not picked up atoms to the left and free spots available
        auto possibleSites = arch.getSitesInRow(zone, row);
        possibleSites.erase(
            std::remove_if(possibleSites.begin(), possibleSites.end(),
                           [&](const auto s) { return !initialFreeSites[s]; }),
            possibleSites.end());
        const auto s =
            possibleSites[std::min(notPickedUpLeft, freeSpotsInRow - 1)];
        placement[q].positionStatus   = Atom::PositionStatus::DEFINED;
        *placement[q].initialPosition = arch.getPositionOfSite(s);
        initialFreeSites[s]           = false;
        currentFreeSites[s]           = false;
      }
      // here the position of q is defined
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      std::vector<std::shared_ptr<Point>> loadStart;
      std::vector<std::shared_ptr<Point>> loadEnd;
      const auto currentX = placement[q].currentPosition->x;
      const auto y        = placement[q].currentPosition->y;
      // pick up q itself
      loadStart.emplace_back(placement[q].currentPosition);
      currentFreeSites[arch.getSiteAt(*placement[q].currentPosition)] = true;
      placement[q].currentPosition = std::make_shared<Point>(currentX + d, y);
      loadEnd.emplace_back(placement[q].currentPosition);
      currentlyShuttling.insert(q);
      // iterate through all not yet picked up atoms to the left and check
      // whether they can be picked up
      auto x = currentX - d;
      if (i > 0) {
        for (std::size_t j = i - 1;; --j) {
          const auto p = moveableOrdered[j];
          // if the atom is picked up,
          // move it to the correct row
          if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
            start.emplace_back(placement[p].currentPosition);
            placement[p].currentPosition = std::make_shared<Point>(x, y);
            end.emplace_back(placement[p].currentPosition);
            const auto nx = arch.getNearestXLeft(x);
            x             = nx == x ? x - dx : nx - d;
          } else {
            // check whether j can be
            // picked up
            if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
              if (placement[p].currentPosition->y == y and
                  placement[p].currentPosition->x <= x - d) {
                // pick up p
                pickUpOrderMoveable.erase(
                    std::remove(pickUpOrderMoveable.begin(),
                                pickUpOrderMoveable.end(), p),
                    pickUpOrderMoveable.end());
                x = placement[p].currentPosition->x;
                loadStart.emplace_back(placement[p].currentPosition);
                currentFreeSites[arch.getSiteAt(
                    *placement[p].currentPosition)] = true;
                placement[p].currentPosition =
                    std::make_shared<Point>(x + d, y);
                loadEnd.emplace_back(placement[p].currentPosition);
                currentlyShuttling.insert(p);
                x -= d;
              }
            } else {
              // if atom is not placed
              // yet find a site to the
              // left that it can be
              // picked up together with
              // q, i.e. find next free
              // site to the left
              std::int64_t freeX = x;
              bool         free  = false;
              try {
                while (!free) {
                  const auto& site = arch.getNearestSiteLeft({freeX, y});
                  freeX            = arch.getPositionOfSite(site).x;
                  if (initialFreeSites[site] and
                      std::find(placement[p].zones.cbegin(),
                                placement[p].zones.cend(),
                                arch.getZoneOfSite(site)) !=
                          placement[p].zones.cend()) {
                    // the site is free
                    // and satisfies the
                    // zone restrictions
                    // of the atom
                    free = true;
                  } else {
                    freeX -= d;
                  }
                }
              } catch (std::invalid_argument& e) {
                // if x reached the left
                // end and there was no
                // site to the left
                // anymore, free remains
                // false in this case
              }
              if (free) {
                // place p on the free
                // site
                placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
                *placement[p].initialPosition = {freeX, y};
                initialFreeSites[arch.getSiteAt(
                    *placement[p].initialPosition)] = false;
                // pick up p
                pickUpOrderMoveable.erase(
                    std::remove(pickUpOrderMoveable.begin(),
                                pickUpOrderMoveable.end(), p),
                    pickUpOrderMoveable.end());
                loadStart.emplace_back(placement[p].currentPosition);
                placement[p].currentPosition =
                    std::make_shared<Point>(freeX + d, y);
                loadEnd.emplace_back(placement[p].currentPosition);
                currentlyShuttling.insert(p);
                x = freeX - d;
              }
            }
          }
          if (j == 0) {
            break;
          }
        }
      }
      // iterate through all not yet picked up atoms to the right and check
      // whether they can be picked up
      x = arch.getNearestXRight(currentX) + d;
      for (std::size_t j = i + 1; j < moveableOrdered.size(); ++j) {
        const auto p = moveableOrdered[j];
        // if the atom is picked up, move it to the correct row
        if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
          start.emplace_back(placement[p].currentPosition);
          placement[p].currentPosition = std::make_shared<Point>(x, y);
          end.emplace_back(placement[p].currentPosition);
          const auto nx = arch.getNearestXRight(x + d);
          x             = nx == x + d ? nx - d + dx : nx + d;
        } else {
          // check whether j can be
          // picked up
          if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
            if (placement[p].currentPosition->y == y and
                placement[p].currentPosition->x >= x - d) {
              // pick up p
              pickUpOrderMoveable.erase(std::remove(pickUpOrderMoveable.begin(),
                                                    pickUpOrderMoveable.end(),
                                                    p),
                                        pickUpOrderMoveable.end());
              x = placement[p].currentPosition->x;
              loadStart.emplace_back(placement[p].currentPosition);
              currentFreeSites[arch.getSiteAt(*placement[p].currentPosition)] =
                  true;
              placement[p].currentPosition = std::make_shared<Point>(x + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
              currentlyShuttling.insert(p);
              const auto nx = arch.getNearestXRight(x + d);
              x             = nx == x + d ? nx - d + dx : nx + d;
            }
          } else {
            // if atom is not placed yet
            // find a site to the right
            // that it can be picked up
            // together with q, i.e. find
            // next free site to the
            // right
            std::int64_t freeX = x;
            bool         free  = false;
            try {
              while (!free) {
                const auto& site = arch.getNearestSiteRight({freeX - d, y});
                freeX            = arch.getPositionOfSite(site).x;
                if (initialFreeSites[site] and
                    std::find(placement[p].zones.cbegin(),
                              placement[p].zones.cend(),
                              arch.getZoneOfSite(site)) !=
                        placement[p].zones.cend()) {
                  // the site is free and
                  // satisfies the zone
                  // restrictions of the
                  // atom
                  free = true;
                } else {
                  freeX = arch.getNearestXRight(freeX + d) + d;
                }
              }
            } catch (std::invalid_argument& e) {
              // if x reached the right
              // end and there was no
              // site to the right
              // anymore, free remains
              // false in this case
            }
            if (free) {
              // place p on the free site
              placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
              *placement[p].initialPosition = {freeX, y};
              initialFreeSites[arch.getSiteAt(*placement[p].initialPosition)] =
                  false;
              // pick up p
              pickUpOrderMoveable.erase(std::remove(pickUpOrderMoveable.begin(),
                                                    pickUpOrderMoveable.end(),
                                                    p),
                                        pickUpOrderMoveable.end());
              loadStart.emplace_back(placement[p].currentPosition);
              placement[p].currentPosition =
                  std::make_shared<Point>(freeX + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
              currentlyShuttling.insert(p);
              const auto nx = arch.getNearestXRight(freeX + d);
              x             = nx == freeX + d ? nx - d + dx : nx + d;
            }
          }
        }
      }
      if (!start.empty()) {
        mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(LOAD, loadStart, loadEnd);
    }
    // all atoms are picked up in order, move them to the interaction zone and
    // store them there
    std::vector<std::shared_ptr<Point>> startMoveable;
    std::vector<std::shared_ptr<Point>> endMoveable;
    std::transform(moveableOrdered.cbegin(), moveableOrdered.cend(),
                   std::back_inserter(startMoveable),
                   [&](const auto& q) { return placement[q].currentPosition; });
    // ------------------------------------------------------------------
    // 4. Apply the cz gates

    for (const auto& timeframe : moveable) {
      for (const auto q : moveableOrdered) {
        const auto x = timeframe.at(q);
        if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
          std::stringstream ss;
          ss << "Atom " << q
             << " was unexpectedly not "
                "picked up.";
          throw std::logic_error(ss.str());
        }
        if (x >= 0 and static_cast<std::size_t>(x) < sites.size()) {
          const auto pos =
              arch.getPositionOfSite(sites[static_cast<std::size_t>(x)]);
          placement[q].currentPosition =
              std::make_shared<Point>(pos.x, pos.y + d);
        } else if (x < 0) {
          const auto pos = arch.getPositionOfSite(sites[0]);
          placement[q].currentPosition =
              std::make_shared<Point>(pos.x + (x * dx), pos.y + d);
        } else { // x >= sites.size()
          const auto pos = arch.getPositionOfSite(sites[sites.size() - 1]);
          placement[q].currentPosition =
              std::make_shared<Point>(pos.x + (x * dx), pos.y + d);
        }
        endMoveable.emplace_back(placement[q].currentPosition);
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, startMoveable,
                                                 endMoveable);
      mappedQc.emplaceBack<NAGlobalOperation>(OpType{qc::OpType::Z, 1});
      for (const auto& q : moveableOrdered) {
        for (const auto& p : fixedOrdered) {
          const auto qPos = *placement[q].currentPosition;
          const auto pPos = *placement[p].currentPosition;
          if ((qPos - pPos).length() <= arch.getInteractionRadius()) {
            graph.getEdge(p, q)->execute();
          }
        }
      }
      startMoveable.clear();
      startMoveable = endMoveable;
      endMoveable.clear();
    }
    // -----------------------------------------------------------------
    // 5. move the atoms back to the storage zone
    std::vector<std::tuple<Zone, Index, std::size_t>> freeSitesPerRow;
    for (const auto& z : initialZones) {
      for (std::size_t r = 0; r < arch.getNrowsInZone(z); ++r) {
        const auto& sitesInRow = arch.getSitesInRow(z, r);
        freeSitesPerRow.emplace_back(
            z, r,
            std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                            [&](const auto& acc, const auto& s) {
                              return acc + (initialFreeSites[s] ? 1 : 0);
                            }));
      }
    }
    std::sort(freeSitesPerRow.begin(), freeSitesPerRow.end(),
              [&](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::size_t moveableSpotsNeeded = moveableOrdered.size();
    std::vector<std::tuple<Zone, Index, std::size_t>> moveableSelectedRows;
    for (std::size_t i = 0; i < freeSitesPerRow.size(); ++i) {
      const auto& [z, r, n] = freeSitesPerRow[i];
      if (n >= moveableSpotsNeeded) {
        if (i + 1 == freeSitesPerRow.size() or
            std::get<2>(freeSitesPerRow[i + 1]) < moveableSpotsNeeded) {
          moveableSelectedRows.emplace_back(z, r, moveableSpotsNeeded);
          moveableSpotsNeeded = 0;
          freeSitesPerRow[i]  = {z, r, n - moveableSpotsNeeded};
          break;
        }
      } else {
        moveableSelectedRows.emplace_back(z, r, n);
        moveableSpotsNeeded -= n;
        freeSitesPerRow[i] = {z, r, 0};
      }
    }
    std::sort(moveableSelectedRows.begin(), moveableSelectedRows.end(),
              [&](const auto& a, const auto& b) {
                return (std::get<0>(a) == std::get<0>(b) and
                        std::get<1>(a) < std::get<1>(b)) or
                       std::get<0>(a) < std::get<0>(b);
              });
    for (auto& [z, r, n] : moveableSelectedRows) {
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      std::vector<std::shared_ptr<Point>> storeStart;
      std::vector<std::shared_ptr<Point>> storeEnd;
      for (std::size_t i = 0; i < moveableOrdered.size(); ++i) {
        const auto& q = moveableOrdered[i];
        start.emplace_back(placement[q].currentPosition);
        if (n > 0) {
          const auto& sitesInRow = arch.getSitesInRow(z, r);
          // find next (leftmost) free
          // site in the row
          const auto& s =
              *std::find_if(sitesInRow.cbegin(), sitesInRow.cend(),
                            [&](const auto& s) { return currentFreeSites[s]; });
          const auto& sPos = arch.getPositionOfSite(s);
          // store the current atom in
          // the free spot either if
          // there are as many atoms as
          // free spots to the right (all
          // can be stored in this line)
          // or the atom is closer than
          // the next one to the right
          if (n == moveableOrdered.size() - i or
              std::abs(placement[q].currentPosition->x - sPos.x) <
                  std::abs(
                      placement[moveableOrdered[i + 1]].currentPosition->x -
                      sPos.x)) {
            placement[q].currentPosition =
                std::make_shared<Point>(sPos.x + d, sPos.y);
            end.emplace_back(placement[q].currentPosition);
            storeStart.emplace_back(placement[q].currentPosition);
            placement[q].currentPosition =
                std::make_shared<Point>(sPos.x, sPos.y);
            storeEnd.emplace_back(placement[q].currentPosition);
            currentlyShuttling.erase(q);
            currentFreeSites[s] = false;
            initialFreeSites[s] = false;
            n -= 1;
          } else {
            placement[q].currentPosition = std::make_shared<Point>(
                arch.getNearestXLeft(placement[q].currentPosition->x) + d,
                sPos.y);
            end.emplace_back(placement[q].currentPosition);
          }
        } else {
          const auto y = arch.getPositionOfSite(arch.getSitesInRow(z, r)[0]).y;
          placement[q].currentPosition = std::make_shared<Point>(
              arch.getNearestXLeft(placement[q].currentPosition->x) + d, y);
          end.emplace_back(placement[q].currentPosition);
        }
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
      mappedQc.emplaceBack<NAShuttlingOperation>(STORE, storeStart, storeEnd);
    }
    // -------------------------------------------------------------
    startFixed.clear();
    endFixed.clear();
    for (const auto& q : fixedOrdered) {
      startFixed.emplace_back(placement[q].currentPosition);
      const auto pos = *placement[q].currentPosition;
      currentlyShuttling.insert(q);
      currentFreeSites[arch.getSiteAt(*placement[q].currentPosition)] = true;
      placement[q].currentPosition = std::make_shared<Point>(pos.x + d, pos.y);
      endFixed.emplace_back(placement[q].currentPosition);
    }
    mappedQc.emplaceBack<NAShuttlingOperation>(LOAD, startFixed, endFixed);
    std::sort(freeSitesPerRow.begin(), freeSitesPerRow.end(),
              [&](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::size_t fixedSpotsNeeded = fixedOrdered.size();
    std::vector<std::tuple<Zone, Index, std::size_t>> fixedSelectedRows;
    for (std::size_t i = 0; i < freeSitesPerRow.size(); ++i) {
      const auto& [z, r, n] = freeSitesPerRow[i];
      if (n >= fixedSpotsNeeded) {
        if (i + 1 == freeSitesPerRow.size() or
            std::get<2>(freeSitesPerRow[i + 1]) < fixedSpotsNeeded) {
          fixedSelectedRows.emplace_back(z, r, fixedSpotsNeeded);
          fixedSpotsNeeded   = 0;
          freeSitesPerRow[i] = {z, r, n - fixedSpotsNeeded};
          break;
        }
      } else {
        fixedSelectedRows.emplace_back(z, r, n);
        fixedSpotsNeeded -= n;
        freeSitesPerRow[i] = {z, r, 0};
      }
    }
    std::sort(fixedSelectedRows.begin(), fixedSelectedRows.end(),
              [&](const auto& a, const auto& b) {
                return (std::get<0>(a) == std::get<0>(b) and
                        std::get<1>(a) < std::get<1>(b)) or
                       std::get<0>(a) < std::get<0>(b);
              });
    for (auto& [z, r, n] : fixedSelectedRows) {
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      std::vector<std::shared_ptr<Point>> storeStart;
      std::vector<std::shared_ptr<Point>> storeEnd;
      for (std::size_t i = 0; i < fixedOrdered.size(); ++i) {
        const auto& q = fixedOrdered[i];
        start.emplace_back(placement[q].currentPosition);
        if (n > 0) {
          const auto& sitesInRow = arch.getSitesInRow(z, r);
          // find next (leftmost) free
          // site in the row
          const auto& s =
              *std::find_if(sitesInRow.cbegin(), sitesInRow.cend(),
                            [&](const auto& s) { return currentFreeSites[s]; });
          const auto& sPos = arch.getPositionOfSite(s);
          // store the current atom in
          // the free spot either if
          // there are as many atoms as
          // free spots to the right (all
          // can be stored in this line)
          // or the atom is closer than
          // the next one to the right
          if (n == fixedOrdered.size() - i or
              std::abs(placement[q].currentPosition->x - sPos.x) <
                  std::abs(placement[fixedOrdered[i + 1]].currentPosition->x -
                           sPos.x)) {
            placement[q].currentPosition =
                std::make_shared<Point>(sPos.x + d, sPos.y);
            end.emplace_back(placement[q].currentPosition);
            storeStart.emplace_back(placement[q].currentPosition);
            placement[q].currentPosition =
                std::make_shared<Point>(sPos.x, sPos.y);
            storeEnd.emplace_back(placement[q].currentPosition);
            currentlyShuttling.erase(q);
            currentFreeSites[s] = false;
            initialFreeSites[s] = false;
            n -= 1;
          } else {
            placement[q].currentPosition = std::make_shared<Point>(
                arch.getNearestXLeft(placement[q].currentPosition->x) + d,
                sPos.y);
            end.emplace_back(placement[q].currentPosition);
          }
        } else {
          const auto y = arch.getPositionOfSite(arch.getSitesInRow(z, r)[0]).y;
          placement[q].currentPosition = std::make_shared<Point>(
              arch.getNearestXLeft(placement[q].currentPosition->x) + d, y);
          end.emplace_back(placement[q].currentPosition);
        }
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
      mappedQc.emplaceBack<NAShuttlingOperation>(STORE, storeStart, storeEnd);
    }
    // -------------------------------------------------------------
    it = (*executableSet)->begin();
  }
  for (auto& p : placement) {
    if (p.positionStatus == Atom::PositionStatus::UNDEFINED) {
      // find next free site in first zone
      std::vector<Index> possibleSites;
      for (const auto& z : p.zones) {
        for (const auto& s : arch.getSitesInZone(z)) {
          possibleSites.emplace_back(s);
        }
      }
      const auto& freeSite =
          std::find_if(possibleSites.cbegin(), possibleSites.cend(),
                       [&](const auto& s) { return initialFreeSites[s]; });
      *p.initialPosition          = arch.getPositionOfSite(*freeSite);
      p.positionStatus            = Atom::PositionStatus::DEFINED;
      initialFreeSites[*freeSite] = false;
      currentFreeSites[*freeSite] = false;
    }
    mappedQc.emplaceInitialPosition(p.initialPosition);
  }
  //========================= END MAPPING =========================
  // get end time
  auto startPostprocess = std::chrono::high_resolution_clock::now();
  postprocess();
  auto end = std::chrono::high_resolution_clock::now();
  // build remaining statistics
  stats.numInitialGates = qc.getNops();
  stats.numMappedGates  = qc.getNops();
  // get the mapping time in milliseconds
  stats.preprocessTime =
      std::chrono::duration<qc::fp, std::milli>(startMapping - startPreprocess).count();
  stats.mappingTime =
      std::chrono::duration<qc::fp, std::milli>(startPostprocess - startMapping).count();
  stats.postprocessTime =
      std::chrono::duration<qc::fp, std::milli>(end - startPostprocess).count();
  done = true;
}
} // namespace na