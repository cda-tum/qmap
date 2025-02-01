//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Definitions.hpp"
#include "datastructures/DirectedAcyclicGraph.hpp"
#include "datastructures/Layer.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

using Color = std::uint16_t;
using Edge = std::pair<qc::Qubit, qc::Qubit>;
using InteractionGraph = qc::Layer::InteractionGraph;

class NAGraphAlgorithms {
public:
  /**
   * @brief Computes some maximal independent set of vertices.
   * @details The method iterates over the nodes in the interaction graph
   * ordered by their degree in descending order. For each node, it checks
   * whether it is adjacent to any node in the current maximal independent set.
   * If not, the node is added to the set. At the end, the set cannot be
   * increased in size by adding any other node without violating the
   * independence property.
   * @param g the graph
   * @return a maximal independent set of vertices
   */
  [[nodiscard]] static auto getMaxIndependentSet(const InteractionGraph& g)
      -> std::unordered_set<qc::Qubit>;

  /**
   * @brief This method returns all edges that are incident to some vertex in
   * the given set of vertices.
   * @param g the graph
   * @param vs the set of vertices
   * @return a set of edges with at least one node in vs
   */
  [[nodiscard]] static auto
  coveredEdges(const InteractionGraph& g,
               const std::unordered_set<qc::Qubit>& vs)
      -> std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>>;

  /**
   * @brief Get the Least Admissible Color for an edge.
   * @details For a coloring to be valid no two adjacent edges can have the same
   * color. Consequently, the least admissible color is one that is not used by
   * any adjacent edge.
   * Additionally, the color must be greater than the maximum color of any
   * adjacent edge that does not contain the vertex v to ensure that the
   * following constraint is satisfied:
   * Let two nodes u and v be in the minimum maximal independent set and the
   * node w and w' both adjacent to u and v. The edge (u, w) has a smaller
   * coloring than the edge (w, v) iff the edge (u, w') has a smaller coloring
   * than the edge (w', v), e.g.
   *
   *                (u)—– 0 —(w)
   *                  \        \
   *                   3        1
   *                    \        \
   *                    (w')— 4 —(v)
   *
   * A '2' instead of the '4' would not be allowed the other colors unchanged.
   *
   * @param coloring a mapping from edges to colors
   * @param maxColor the maximum color used so far
   * @param e the edge to be colored
   * @param v the root of the edge, i.e. the endpoint of the edge that is
   * contained in the previously computed maximal independent set
   * @param sequence the sequence of vertices selected to be the moveable ones
   * @param partialOrder the partial order of the fixed vertices induced by the
   * coloring so far
   * @param ranks the ranks of the vertices in the sequence, i.e., a mapping of
   * nodes in the sequence to their index in the sequence, i.e., their rank
   * @return the least admissible color
   */
  [[nodiscard]] static auto getLeastAdmissibleColor(
      const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
          coloring,
      const Color& maxColor, const Edge& e, const qc::Qubit& v,
      const std::vector<qc::Qubit>& sequence,
      const qc::DirectedAcyclicGraph<qc::Qubit>& partialOrder,
      const std::unordered_map<std::pair<qc::Qubit, Color>, std::size_t,
                               qc::PairHash<qc::Qubit, Color>>& ranks) -> Color;

  /**
   * @brief Colors all given edges starting with edges that are adjacent to the
   * first vertex in the queue.
   * @param g the interaction graph
   * @param edges e.g. edges adjacent to a selected vertex
   * @param nodesQueue could for example be sorted by degree, highest first
   * @return a valid coloring of the edges
   */
  [[nodiscard]] static auto colorEdges(
      const InteractionGraph& g,
      const std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>>& edges,
      const std::vector<qc::Qubit>& nodesQueue)
      -> std::pair<
          std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>,
          qc::DirectedAcyclicGraph<qc::Qubit>>;

  /**
   * @brief Computes the resting positions of the moveable vertices between the
   * fixed vertices.
   * @param moveable all moveable vertices in order
   * @param fixed all fixed vertices in order
   * @param coloring the coloring determining the interaction of each moveable
   * and fixed vertex in every time step
   * @return the required places for the moveable vertices to rest in time steps
   * when they do not undergo an interaction
   */
  [[nodiscard]] static auto computeRestingPositions(
      const std::vector<qc::Qubit>& moveable,
      const std::vector<qc::Qubit>& fixed,
      const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
          coloring) -> std::vector<std::size_t>;

  /**
   * @brief Groups the vertices in the sequence by connected components.
   * @param g the graph with connected components
   * @param sequence the sequence of vertices
   * @return a new sequence with the vertices grouped by connected components
   */
  [[nodiscard]] static auto
  groupByConnectedComponent(const InteractionGraph& g,
                            const std::vector<qc::Qubit>& sequence)
      -> std::vector<qc::Qubit>;

  /**
   * @brief The result will be (1) a mapping of vertices to their relative x
   * position in every time step (the moveable vertices) and (2) a mapping of
   * vertices to their x position (at all time step, the fixed vertices).
   * @details The algorithm works as follows:
   * 1. Compute a maximal independent set of vertices. This partitions the
   * vertices into moveable and fixed vertices.
   * 2. Color the edges of the interaction graph such that no two adjacent edges
   * have the same color.
   * 3. Find relative x-positions for the fixed vertices.
   * 4. Compute the resting positions of the moveable vertices between the fixed
   * vertices.
   * @param g        the interaction graph
   * @param maxSites the number of sites horizontally in the entnagling zones
   * @return A tuple of vectors containing the moveable vertices, the fixed
   * vertices mapped to their relative x-position in every time step.
   */
  [[nodiscard]] static auto computeSequence(const InteractionGraph& g,
                                            std::size_t maxSites)
      -> std::pair<std::vector<std::unordered_map<qc::Qubit, std::int64_t>>,
                   std::unordered_map<qc::Qubit, std::int64_t>>;
};
} // namespace na
