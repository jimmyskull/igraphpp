#ifndef IGRAPHPP_IGRAPH_HPP_
#define IGRAPHPP_IGRAPH_HPP_

#include <utility>

#include <igraph/igraph.h>

#include "./exception.hpp"

namespace igraph {

using Edge = std::pair<int, int>;

static constexpr double kInfinity = IGRAPH_INFINITY;

static constexpr bool Directed = IGRAPH_DIRECTED;
static constexpr bool Undirected = IGRAPH_UNDIRECTED;

enum class Connectedness { Weak = IGRAPH_WEAK, Strong = IGRAPH_STRONG };

enum class MultiEdges { Ignore = false, OneOf = true };

enum class Mode { Out = IGRAPH_OUT, In = IGRAPH_IN, All = IGRAPH_ALL };

enum class Mutuality { None = false, MutualConnections = true };

enum class Periodicity { None = false, Circular = true };

enum class Loops { None = IGRAPH_NO_LOOPS, Allow = IGRAPH_LOOPS };

enum class EdgeOrder {
  ById = IGRAPH_EDGEORDER_ID,
  BySourceId = IGRAPH_EDGEORDER_FROM,
  ByTargetId = IGRAPH_EDGEORDER_TO
};

enum class AdjacencyMatrixMode {
  Directed = IGRAPH_ADJ_DIRECTED,
  Undirected = IGRAPH_ADJ_UNDIRECTED,
  Max = IGRAPH_ADJ_MAX,
  Min = IGRAPH_ADJ_MIN,
  Plus = IGRAPH_ADJ_PLUS,
  Upper = IGRAPH_ADJ_UPPER,
  Lower = IGRAPH_ADJ_LOWER
};

enum class StarMode {
  Out = IGRAPH_STAR_OUT,
  In = IGRAPH_STAR_IN,
  Mutual = IGRAPH_STAR_MUTUAL,
  Undirected = IGRAPH_STAR_UNDIRECTED
};

enum class TreeMode {
  Out = IGRAPH_TREE_OUT,
  In = IGRAPH_TREE_IN,
  Undirected = IGRAPH_TREE_UNDIRECTED
};

enum class BarabasiAlgorithm {
  Bag = IGRAPH_BARABASI_BAG,
  PSumTree = IGRAPH_BARABASI_PSUMTREE,
  PSumTreeMultiple = IGRAPH_BARABASI_PSUMTREE_MULTIPLE
};

enum class DegreeSequenceMethod {
  Simple = IGRAPH_DEGSEQ_SIMPLE,
  SimpleNoMultiple = IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE,
  VL = IGRAPH_DEGSEQ_VL
};

enum class RewiringMode {
  Simple = IGRAPH_REWIRING_SIMPLE,
  SimpleLoops = IGRAPH_REWIRING_SIMPLE_LOOPS
};

enum class SubgraphImplementation {
  CopyAndDelete = IGRAPH_SUBGRAPH_COPY_AND_DELETE,
  CreateFromScratch = IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH
};

enum class DirectedMode {
  Arbitrary = IGRAPH_TO_DIRECTED_ARBITRARY,
  Mutual = IGRAPH_TO_DIRECTED_MUTUAL
};

enum class UndirectedMode {
  Each = IGRAPH_TO_UNDIRECTED_EACH,
  Collapse = IGRAPH_TO_UNDIRECTED_COLLAPSE,
  Mutual = IGRAPH_TO_UNDIRECTED_MUTUAL
};

int SetSeed(unsigned long int seed);

}  // namespace igraph

#include "./mapper.hpp"
#include "./vectorptr.hpp"
#include "./vector.hpp"
#include "./vertex_selector.hpp"
#include "./vertex_iterator.hpp"
#include "./edge_selector.hpp"
#include "./edge_iterator.hpp"
#include "./matrix.hpp"
#include "./graph.hpp"

#endif  // IGRAPHPP_IGRAPH_HPP_
