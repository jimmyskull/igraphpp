#ifndef IGRAPHPP_IGRAPH_HPP_
#define IGRAPHPP_IGRAPH_HPP_

#include <utility>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

using Edge = std::pair<int, int>;

static constexpr double kInfinity = IGRAPH_INFINITY;

enum Directedness {
  Undirected = IGRAPH_UNDIRECTED,
  Directed = IGRAPH_DIRECTED
};

enum Connectedness {
  WeaklyConnected = IGRAPH_WEAK,
  StronglyConnected = IGRAPH_STRONG
};

enum MultiEdges { IgnoreMultiEdges = false, OneOfMultiEdges = true };

enum NeighborMode { Out = IGRAPH_OUT, In = IGRAPH_IN, All = IGRAPH_ALL };

enum Mutuality { NotMutual = false, MutualConnections = true };

enum Periodicity { NotPeriodic = false, Circular = true };

enum Loops { NoLoops = IGRAPH_NO_LOOPS, AllowLoops = IGRAPH_LOOPS };

enum EdgeOrder {
  EdgeById = IGRAPH_EDGEORDER_ID,
  EdgeBySourceId = IGRAPH_EDGEORDER_FROM,
  EdgeByTargetId = IGRAPH_EDGEORDER_TO
};

enum AdjacencyMatrixMode {
  AdjacencyDirected = IGRAPH_ADJ_DIRECTED,
  AdjacencyUndirected = IGRAPH_ADJ_UNDIRECTED,
  AdjacencyMax = IGRAPH_ADJ_MAX,
  AdjacencyMin = IGRAPH_ADJ_MIN,
  AdjacencyPlus = IGRAPH_ADJ_PLUS,
  AdjacencyUpper = IGRAPH_ADJ_UPPER,
  AdjacencyLower = IGRAPH_ADJ_LOWER
};

enum StarMode {
  StarOut = IGRAPH_STAR_OUT,
  StarIn = IGRAPH_STAR_IN,
  StarMutual = IGRAPH_STAR_MUTUAL,
  StarUndirected = IGRAPH_STAR_UNDIRECTED
};

enum TreeMode {
  TreeOut = IGRAPH_TREE_OUT,
  TreeIn = IGRAPH_TREE_IN,
  TreeUndirected = IGRAPH_TREE_UNDIRECTED
};

enum BarabasiAlgorithm {
  BarabasiBag = IGRAPH_BARABASI_BAG,
  BarabasiPSumTree = IGRAPH_BARABASI_PSUMTREE,
  BarabasiPSumTreeMultiple = IGRAPH_BARABASI_PSUMTREE_MULTIPLE
};

enum DegreeSequenceMethod {
  DegreeSequenceSimple = IGRAPH_DEGSEQ_SIMPLE,
  DegreeSequenceSimpleNoMultiple = IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE,
  DegreeSequenceVL = IGRAPH_DEGSEQ_VL
};

enum RewiringMode {
  RewiringSimple = IGRAPH_REWIRING_SIMPLE,
  RewiringSimpleLoops = IGRAPH_REWIRING_SIMPLE_LOOPS
};

enum SubgraphImplementation {
  SubgraphCopyAndDelete = IGRAPH_SUBGRAPH_COPY_AND_DELETE,
  SubgraphCreateFromScratch = IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH
};

enum DirectedMode {
  DirectedArbitrary = IGRAPH_TO_DIRECTED_ARBITRARY,
  DirectedMutual = IGRAPH_TO_DIRECTED_MUTUAL
};

enum UndirectedMode {
  UndirectedEach = IGRAPH_TO_UNDIRECTED_EACH,
  UndirectedCollapse = IGRAPH_TO_UNDIRECTED_COLLAPSE,
  UndirectedMutual = IGRAPH_TO_UNDIRECTED_MUTUAL
};

static inline int SetSeed(unsigned long int seed) {
  igraph_rng_t *rng = igraph_rng_default();
  int ret = igraph_rng_seed(rng, seed);
  return SafeCall(ret);
}

} // namespace igraph

#include "./mapper.hpp"
#include "./vector.hpp"
#include "./vectorptr.hpp"
#include "./vertex_selector.hpp"
#include "./vertex_iterator.hpp"
#include "./edge_selector.hpp"
#include "./edge_iterator.hpp"
#include "./matrix.hpp"
#include "./graph.hpp"

#endif // IGRAPHPP_IGRAPH_HPP_
