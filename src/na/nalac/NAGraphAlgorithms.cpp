#include "na/nalac/NAGraphAlgorithms.hpp"

#include "Definitions.hpp"
#include "datastructures/DirectedAcyclicGraph.hpp"
#include "datastructures/DisjointSet.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <list>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

auto NAGraphAlgorithms::getMaxIndependentSet(const InteractionGraph& g)
    -> std::unordered_set<qc::Qubit> {
  std::unordered_set<qc::Qubit> result;
  const auto& vertices = g.getVertices();
  std::list<qc::Qubit> queue(vertices.cbegin(), vertices.cend());
  // sort the vertices by degree in descending order
  queue.sort([&](const auto& u, const auto& v) {
    return g.getDegree(u) > g.getDegree(v);
  });
  while (!queue.empty()) {
    const qc::Qubit v = queue.front();
    result.emplace(v);
    queue.remove_if(
        [&](const qc::Qubit& u) { return u == v or g.isAdjacent(u, v); });
  }
  return result;
}

auto NAGraphAlgorithms::coveredEdges(const InteractionGraph& g,
                                     const std::unordered_set<qc::Qubit>& vs)
    -> std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>> {
  std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>> result;
  for (const qc::Qubit& v : vs) {
    try {
      result.merge(g.getAdjacentEdges(v));
    } catch (const std::invalid_argument&) {
      std::stringstream ss;
      ss << "Vertices must be a subset of all vertices in the graph, in "
            "particular "
         << v << " is not in the graph.";
      throw std::invalid_argument(ss.str());
    }
  }
  return result;
}

auto NAGraphAlgorithms::getLeastAdmissibleColor(
    const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
        coloring,
    const Color& maxColor, const Edge& e, const qc::Qubit& v,
    const std::vector<qc::Qubit>& sequence,
    const qc::DirectedAcyclicGraph<qc::Qubit>& partialOrder,
    const std::unordered_map<std::pair<qc::Qubit, Color>, std::size_t,
                             qc::PairHash<qc::Qubit, Color>>& ranks) -> Color {
  // compute the minimum admissible color as the maximum color +1 of adjacent
  // edges that do not contain the vertex v
  Color minAdmissibleColor = 0;
  for (const auto& [f, k] : coloring) {
    if (InteractionGraph::isAdjacentEdge(e, f) and v != f.first and
        v != f.second) {
      minAdmissibleColor =
          std::max(minAdmissibleColor, static_cast<Color>(k + 1));
    }
  }

  std::vector freeColors(maxColor + 2U - minAdmissibleColor, true);
  for (const auto& [f, k] : coloring) {
    if (InteractionGraph::isAdjacentEdge(e, f) && k >= minAdmissibleColor &&
        k <= maxColor + 1) {
      freeColors[k - minAdmissibleColor] = false;
    }
  }
  // The minAdmissibleColor is now the minimum of the free colors
  // that does not generate a cycle in the graph induced by the partial order
  // of SLM traps.

  // u is the SLM trap that is adjacent to v in the edge e
  const qc::Qubit u = v == e.first ? e.second : e.first;
  // get the rank of u, i.e., the index of v in the sequence of AOD traps
  // (not contained in ranks since the edge e is not colored yet)
  const auto rankOfU = static_cast<std::size_t>(std::distance(
      sequence.cbegin(), std::find(sequence.cbegin(), sequence.cend(), v)));
  for (std::size_t i = 0; i < freeColors.size(); ++i) {
    if (!freeColors[i]) {
      continue;
    }
    const auto leastAdmissibleColor =
        static_cast<Color>(i + minAdmissibleColor);
    bool isAdmissible = true;
    for (const auto& [f, k] : coloring) {
      if (f.first == v or f.second == v) {
        const qc::Qubit w = f.first == v ? f.second : f.first;
        if (k > leastAdmissibleColor) {
          if (partialOrder.isReachable(w, u)) {
            isAdmissible = false;
            break;
          }
        } else if (k < leastAdmissibleColor) {
          if (partialOrder.isReachable(u, w)) {
            throw std::logic_error("Coloring cannot be completed to a valid "
                                   "one (cycle is unavoidable).");
          }
        }
      } else if (k == leastAdmissibleColor) {
        // get the SLM atom from the edge f
        const qc::Qubit w = std::find(sequence.cbegin(), sequence.cend(),
                                      f.first) == sequence.end()
                                ? f.first
                                : f.second;

        if ((rankOfU > ranks.at({w, k}) && partialOrder.isReachable(w, u)) ||
            (rankOfU < ranks.at({w, k}) && partialOrder.isReachable(u, w))) {
          isAdmissible = false;
          break;
        }
      }
    }
    if (isAdmissible) {
      return leastAdmissibleColor;
    }
  }
  throw std::logic_error("No admissible color found (should never occur).");
}

auto NAGraphAlgorithms::colorEdges(
    const InteractionGraph& g,
    const std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>>& edges,
    const std::vector<qc::Qubit>& nodesQueue)
    -> std::pair<
        std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>,
        qc::DirectedAcyclicGraph<qc::Qubit>> {
  std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>
      coloring{};
  Color maxColor = 0;
  // number of distinct colors of edges adjacent to an edge
  std::unordered_map<Edge, std::size_t, qc::PairHash<qc::Qubit, qc::Qubit>>
      nAdjColors;
  // the degree of the edge seen as a node
  std::unordered_map<Edge, std::size_t, qc::PairHash<qc::Qubit, qc::Qubit>>
      edgeDegree;
  for (const auto& e : edges) {
    nAdjColors[e] = 0;
    edgeDegree[e] = static_cast<std::size_t>(
        std::count_if(edges.cbegin(), edges.cend(), [&](const Edge& f) {
          return InteractionGraph::isAdjacentEdge(e, f);
        }));
  }
  // represent the partial order on the SLM atoms as a directed acyclic graph
  qc::DirectedAcyclicGraph<qc::Qubit> partialOrder;
  const std::unordered_set aodTraps(nodesQueue.cbegin(), nodesQueue.cend());
  for (const auto& v : g.getVertices()) {
    if (aodTraps.find(v) == aodTraps.end()) {
      // only add SLM traps, AOD traps are the one in the queue
      partialOrder.addVertex(v);
    }
  }
  // for every SLM trap and color, calculate its rank, i.e., the index of the
  // interaction partner in the sequence of AOD traps
  std::unordered_map<std::pair<qc::Qubit, Color>, std::size_t,
                     qc::PairHash<qc::Qubit, Color>>
      ranks;

  for (const auto& v : nodesQueue) {
    std::vector<Edge> adjacentEdges{};
    std::copy_if(edges.cbegin(), edges.cend(),
                 std::back_inserter(adjacentEdges),
                 [&](const Edge& e) { return e.first == v or e.second == v; });
    std::sort(adjacentEdges.begin(), adjacentEdges.end(),
              [&](const Edge& a, const Edge& b) {
                const auto u = a.first == v ? a.second : a.first;
                const auto w = b.first == v ? b.second : b.first;
                if (u == w) {
                  // the compare function defines a proper less than relation,
                  // i.e., equal elements must return false
                  return false;
                }
                if (partialOrder.isReachable(u, w)) {
                  // if w is reachable from u, then w needs to come after u,
                  // i.e., u < w
                  return true;
                }
                if (partialOrder.isReachable(w, u)) {
                  // if u is reachable from w, then u needs to come after w,
                  // i.e., w < u and not u < w
                  return false;
                }
                // here: neither u is reachable from w nor w is reachable from u
                if (nAdjColors[a] > nAdjColors[b]) {
                  // this favours nodes with higher number of distinct colors
                  // also see the original implementation of DSatur
                  return true;
                }
                if (nAdjColors[a] < nAdjColors[b]) {
                  // this is just the opposite of the above case
                  return false;
                }
                // Here: nAdjColors[a] == nAdjColors[b]
                if (edgeDegree[a] > edgeDegree[b]) {
                  return true;
                }
                if (edgeDegree[a] < edgeDegree[b]) {
                  return false;
                }
                return u < w;
                // the last line together with the first clause (u != w) is
                // necessary for a well-defined compare function to handle edges
                // that compare equally correctly
              });
    for (const auto& e : adjacentEdges) {
      // color the edge
      coloring[e] = getLeastAdmissibleColor(coloring, maxColor, e, v,
                                            nodesQueue, partialOrder, ranks);
      // update partial order
      const qc::Qubit u = e.first == v ? e.second : e.first;
      ranks[{u, coloring[e]}] = static_cast<std::size_t>(
          std::distance(nodesQueue.cbegin(),
                        std::find(nodesQueue.cbegin(), nodesQueue.cend(), v)));
      for (const auto& [f, k] : coloring) {
        if (f.first == v or f.second == v) {
          const qc::Qubit w = f.first == v ? f.second : f.first;
          if (k < coloring[e]) {
            partialOrder.addEdge(w, u);
          } else if (k > coloring[e]) {
            partialOrder.addEdge(u, w);
          }
        } else if (k == coloring[e]) {
          const qc::Qubit w = std::find(nodesQueue.cbegin(), nodesQueue.cend(),
                                        f.first) == nodesQueue.cend()
                                  ? f.first
                                  : f.second;
          if (ranks[{u, k}] < ranks[{w, k}]) {
            partialOrder.addEdge(w, u);
          } else if (ranks[{u, k}] > ranks[{w, k}]) {
            partialOrder.addEdge(u, w);
          } else {
            throw std::logic_error("Coloring is not valid.");
          }
        }
      }

      maxColor = std::max(maxColor, coloring[e]);
      // update the number of distinct colors of adjacent edges
      for (const auto& f : edges) {
        if (InteractionGraph::isAdjacentEdge(e, f)) {
          std::vector<bool> usedColors(maxColor + 1, false);
          for (const auto& h : edges) {
            if (InteractionGraph::isAdjacentEdge(f, h) and
                coloring.find(h) != coloring.end()) {
              usedColors.at(coloring[h]) = true;
            }
          }
          nAdjColors[f] = static_cast<std::size_t>(
              std::count_if(usedColors.cbegin(), usedColors.cend(),
                            [](const bool& b) { return b; }));
        }
      }
    }
  }
  return {coloring, partialOrder};
}

auto NAGraphAlgorithms::computeRestingPositions(
    const std::vector<qc::Qubit>& moveable, const std::vector<qc::Qubit>& fixed,
    const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
        coloring) -> std::vector<std::size_t> {
  const Color maxColor = std::accumulate(
      coloring.cbegin(), coloring.cend(), static_cast<Color>(0),
      [](Color acc, const auto& value) { return std::max(acc, value.second); });
  std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                     qc::PairHash<std::size_t, std::size_t>>
      resting{};
  for (Color t = 0; t <= maxColor; ++t) {
    // resting positions for timestamp t
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                       qc::PairHash<std::size_t, std::size_t>>
        tResting{};
    // x-coordinates of movable vertices at timestamp t
    std::unordered_map<qc::Qubit, std::size_t> moveableXs{};
    for (const qc::Qubit v : moveable) {
      // get the neighbors of v
      const std::set<qc::Qubit> neighborsOfV = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            const Edge e = value.first;
            if (value.second == t) {
              if (e.first == v) {
                acc.emplace(e.second);
              } else if (e.second == v) {
                acc.emplace(e.first);
              }
            }
            return acc;
          });
      if (!neighborsOfV.empty()) {
        assert(neighborsOfV.size() == 1);
        const qc::Qubit u = *neighborsOfV.begin();
        // get index of u in fixed which is the x-coordinate \wo resting
        const auto& it = std::find(fixed.cbegin(), fixed.cend(), u);
        assert(it != fixed.cend());
        // set x-pos of v to x-pos of u, which is the index of u in fixed
        moveableXs[v] =
            static_cast<std::size_t>(std::distance(fixed.cbegin(), it));
      }
    }
    // map the keys of moveableXs to their indices in moveable
    std::vector<std::size_t> moveableXsIds{};
    std::transform(
        moveableXs.cbegin(), moveableXs.cend(),
        std::back_inserter(moveableXsIds), [&](const auto& value) {
          const auto& it =
              std::find(moveable.cbegin(), moveable.cend(), value.first);
          return static_cast<std::size_t>(std::distance(moveable.cbegin(), it));
        });
    for (const qc::Qubit v : moveable) {
      if (moveableXs.find(v) == moveableXs.end()) {
        // index of v in moveable
        const auto& it = std::find(moveable.cbegin(), moveable.cend(), v);
        const auto i =
            static_cast<std::size_t>(std::distance(moveable.cbegin(), it));
        // get the index of the left neighbor in the keys of moveableXs
        std::vector<std::size_t> leftNeighbors;
        std::copy_if(moveableXsIds.cbegin(), moveableXsIds.cend(),
                     std::back_inserter(leftNeighbors),
                     [&](const auto& j) { return j > i; });
        std::vector<std::size_t> rightNeighbors;
        std::copy_if(moveableXsIds.cbegin(), moveableXsIds.cend(),
                     std::back_inserter(rightNeighbors),
                     [&](const auto& j) { return j < i; });
        // if there is a left and right neighbor, we have to add a
        // resting position
        if (!leftNeighbors.empty() and !rightNeighbors.empty()) {
          // get min left and max right neighbor
          const auto& leftNeighbor =
              std::min_element(leftNeighbors.cbegin(), leftNeighbors.cend());
          const auto& rightNeighbor =
              std::max_element(rightNeighbors.cbegin(), rightNeighbors.cend());
          const auto& pair = std::pair<std::int64_t, std::int64_t>(
              moveableXs[moveable.at(*leftNeighbor)],
              moveableXs[moveable.at(*rightNeighbor)]);
          if (tResting.find(pair) == tResting.end()) {
            tResting[pair] = 1;
          } else {
            tResting[pair]++;
          }
        }
      }
    }
    // tResting is completely computed until here
    // merge tResting with resting
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                       qc::PairHash<std::size_t, std::size_t>>
        newResting;
    for (const auto& elem : resting) {
      const auto& pair = elem.first;
      for (std::size_t i = 0; i < elem.second; ++i) {
        // check whether existing resting positions overlap with newly required
        // ones
        const auto& overlappingPairs =
            std::accumulate(tResting.cbegin(), tResting.cend(),
                            std::vector<std::pair<std::size_t, std::size_t>>(),
                            [&](auto& acc, const auto& value) {
                              if (value.first.first < pair.second &&
                                  pair.first < value.first.second) {
                                acc.emplace_back(value.first);
                              }
                              return acc;
                            });
        if (!overlappingPairs.empty()) {
          // select the smallest one to not loose a lot flexibility
          const std::pair<std::size_t, std::size_t>& minSuperPair =
              *(std::min_element(
                  overlappingPairs.cbegin(), overlappingPairs.cend(),
                  [](const auto& a, const auto& b) {
                    return a.second - a.first < b.second - b.first;
                  }));
          --tResting[minSuperPair];
          if (tResting[minSuperPair] == 0) {
            tResting.erase(minSuperPair);
          }
          const std::pair<std::size_t, std::size_t> newPair{
              std::max(pair.first, minSuperPair.first),
              std::min(pair.second, minSuperPair.second)};
          if (newResting.find(newPair) == newResting.end()) {
            newResting[newPair] = 1;
          } else {
            ++newResting[newPair];
          }
        } else {
          if (newResting.find(pair) == newResting.end()) {
            newResting[pair] = 1;
          } else {
            ++newResting[pair];
          }
        }
      }
    }
    for (const auto& [pair, s] : tResting) {
      if (newResting.find(pair) == newResting.end()) {
        newResting[pair] = s;
      } else {
        newResting[pair] += s;
      }
    }
    resting.clear();
    resting = newResting;
  }
  std::vector<std::size_t> restingPositions{};
  for (const auto& [pair, s] : resting) {
    for (std::size_t i = 0; i < s; ++i) {
      restingPositions.emplace_back(pair.first);
    }
  }
  std::sort(restingPositions.begin(), restingPositions.end());
  return restingPositions;
}

auto NAGraphAlgorithms::groupByConnectedComponent(
    const InteractionGraph& g, const std::vector<qc::Qubit>& sequence)
    -> std::vector<qc::Qubit> {
  const auto& vertices = g.getVertices();
  qc::DisjointSet<qc::Qubit> ds(vertices.cbegin(), vertices.cend());
  for (const qc::Qubit& v : vertices) {
    for (const qc::Qubit& u : vertices) {
      if (g.isAdjacent(v, u)) {
        ds.unionSet(v, u);
      }
    }
  }
  std::vector<qc::Qubit> result{};
  for (const qc::Qubit& v : vertices) {
    if (ds.findSet(v) == v) {
      for (const qc::Qubit& u : sequence) {
        if (ds.findSet(u) == v) {
          result.emplace_back(u);
        }
      }
    }
  }
  return result;
}

auto NAGraphAlgorithms::computeSequence(const InteractionGraph& g,
                                        const std::size_t maxSites)
    -> std::pair<std::vector<std::unordered_map<qc::Qubit, std::int64_t>>,
                 std::unordered_map<qc::Qubit, std::int64_t>> {
  const auto& maxIndepSet = getMaxIndependentSet(g);
  std::vector sequenceUngrouped(maxIndepSet.cbegin(), maxIndepSet.cend());
  // sort the vertices by degree in descending order
  std::sort(sequenceUngrouped.begin(), sequenceUngrouped.end(),
            [&](const auto& u, const auto& v) {
              return g.getDegree(u) > g.getDegree(v);
            });
  auto sequence = groupByConnectedComponent(g, sequenceUngrouped);
  const auto& colorEdgesResult =
      colorEdges(g, coveredEdges(g, maxIndepSet), sequence);
  auto coloring = colorEdgesResult.first;
  auto partialOrder = colorEdgesResult.second;
  auto fixed = partialOrder.orderTopologically();
  auto resting = computeRestingPositions(sequence, fixed, coloring);
  // compute relative x positions of fixed vertices
  std::unordered_map<qc::Qubit, std::int64_t> fixedPositions{};
  for (std::uint32_t x = 0, i = 0; x < fixed.size(); ++x) {
    fixedPositions[fixed[x]] = static_cast<std::int64_t>(x) + i;
    while (i < resting.size() and x == resting[i]) {
      ++i;
    }
  }

  const auto maxSiteUsed =
      std::max_element(
          fixedPositions.cbegin(), fixedPositions.cend(),
          [](const auto& a, const auto& b) { return a.second < b.second; })
          ->second;
  const auto maxSitesSigned = static_cast<std::int64_t>(maxSites);
  if (maxSiteUsed >= maxSitesSigned) {
    // Handle the situation when the entangling zone is not big enough to fit
    // all fixed qubits
    for (auto it = fixedPositions.begin(); it != fixedPositions.end();) {
      if (it->second >= maxSitesSigned) {
        it = fixedPositions.erase(it); // erase returns the next iterator
      } else {
        ++it; // move to the next element
      }
    }
    fixed.erase(std::remove_if(fixed.begin(), fixed.end(),
                               [&fixedPositions](const auto& q) {
                                 return fixedPositions.find(q) ==
                                        fixedPositions.end();
                               }),
                fixed.end());
    for (auto it = coloring.begin(); it != coloring.end();) {
      if (fixedPositions.find(it->first.first) == fixedPositions.end() &&
          fixedPositions.find(it->first.second) == fixedPositions.end()) {
        it = coloring.erase(it); // erase returns the next iterator
      } else {
        ++it; // move to the next element
      }
    }
    sequence.erase(std::remove_if(sequence.begin(), sequence.end(),
                                  [&coloring](const auto& q) {
                                    return !std::any_of(
                                        coloring.cbegin(), coloring.cend(),
                                        [q](const auto& elem) {
                                          return elem.first.first == q ||
                                                 elem.first.second == q;
                                        });
                                  }),
                   sequence.end());
    // recalculate resting positions
    resting = computeRestingPositions(sequence, fixed, coloring);
  }

  // compute relative x positions of moveable vertices at every timestamp
  const Color maxColor = std::accumulate(
      coloring.cbegin(), coloring.cend(), static_cast<Color>(0),
      [](Color acc, const auto& value) { return std::max(acc, value.second); });
  std::vector<std::unordered_map<qc::Qubit, std::int64_t>> moveablePositions(
      maxColor + 1);
  for (Color t = 0; t <= maxColor; ++t) {
    for (const qc::Qubit v : sequence) {
      // get the neighbors of v at timestamp t
      const std::set<qc::Qubit>& neighborsOfV = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            if (value.second == t) {
              const Edge e = value.first;
              if (e.first == v) {
                acc.emplace(e.second);
              } else if (e.second == v) {
                acc.emplace(e.first);
              }
            }
            return acc;
          });
      if (!neighborsOfV.empty()) {
        assert(neighborsOfV.size() == 1);
        const qc::Qubit u = *neighborsOfV.begin();
        // get x-coordinate of u which is a fixed vertex
        const auto& it = fixedPositions.find(u);
        assert(it != fixedPositions.end());
        moveablePositions.at(t)[v] = it->second;
      }
    }
    for (std::size_t i = 0; i < sequence.size(); ++i) {
      if (const qc::Qubit v = sequence[i];
          moveablePositions.at(t).find(v) == moveablePositions.at(t).end()) {
        if (i > 0) {
          // right neighbor already has a position
          const auto rightNeighbor = i - 1;
          const auto rightNeighborX =
              moveablePositions.at(t)[sequence.at(rightNeighbor)];
          const std::int64_t minX =
              std::min(rightNeighborX - 1, static_cast<std::int64_t>(-1));
          assert(rightNeighborX - minX >= 0);
          std::vector<std::int64_t> xrange(
              static_cast<std::size_t>(rightNeighborX - minX));
          std::iota(xrange.begin(), xrange.end(), minX);
          std::vector<std::int64_t> freeX;
          std::copy_if(xrange.cbegin(), xrange.cend(),
                       std::back_inserter(freeX), [&](const std::int64_t x) {
                         return std::all_of(fixedPositions.cbegin(),
                                            fixedPositions.cend(),
                                            [&](const auto& value) {
                                              return x != value.second;
                                            });
                       });
          const auto& maxIt = std::max_element(freeX.cbegin(), freeX.cend());
          moveablePositions.at(t)[v] = *maxIt;
        } else {
          // it is the rightmost atom (w/o position)
          const auto leftNeighbor = std::accumulate(
              moveablePositions.at(t).cbegin(), moveablePositions.at(t).cend(),
              std::make_pair(0UL, -1LL),
              [](const std::pair<const qc::Qubit, std::int64_t>& acc,
                 const auto& value) {
                if (value.second > acc.second) {
                  return value;
                }
                return acc;
              });
          // k is the number of spots between the left positioned neighbor and
          // the leftmost atom (including itself)
          const auto k = static_cast<std::size_t>(std::distance(
              sequence.cbegin(), std::find(sequence.cbegin(), sequence.cend(),
                                           leftNeighbor.first)));
          const std::int64_t maxX =
              std::accumulate(fixedPositions.cbegin(), fixedPositions.cend(),
                              static_cast<std::int64_t>(0),
                              [](const auto& acc, const auto& value) {
                                return std::max(acc, value.second);
                              });
          std::vector<std::int64_t> xrange(
              static_cast<std::size_t>(maxX - leftNeighbor.second));
          std::iota(xrange.begin(), xrange.end(), leftNeighbor.second + 1);
          std::vector<std::int64_t> freeX;
          std::copy_if(xrange.cbegin(), xrange.cend(),
                       std::back_inserter(freeX), [&](const std::int64_t x) {
                         return std::all_of(fixedPositions.cbegin(),
                                            fixedPositions.cend(),
                                            [&](const auto& value) {
                                              return x != value.second;
                                            });
                       });
          if (k <= freeX.size()) {
            moveablePositions.at(t)[v] = freeX[k - 1];
          } else {
            moveablePositions.at(t)[v] =
                maxX + static_cast<std::int64_t>(k) -
                static_cast<std::int64_t>(freeX.size());
          }
        }
      }
    }
  }
  return {moveablePositions, fixedPositions};
}
} // namespace na
