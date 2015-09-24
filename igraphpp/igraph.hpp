#ifndef IGRAPHPP_IGRAPH_HPP_
#define IGRAPHPP_IGRAPH_HPP_

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

enum Directedness {
  Undirected = IGRAPH_UNDIRECTED,
  Directed = IGRAPH_DIRECTED
};

enum Connectedness {
  WeaklyConnected = IGRAPH_WEAK,
  StronglyConnected = IGRAPH_STRONG
};

enum NeighborMode { Out = IGRAPH_OUT, In = IGRAPH_IN, All = IGRAPH_ALL };

enum Mutuality { NotMutual = false, MutualConnections = true };

enum Periodicity { NotPeriodic = false, Circular = true };

enum Loops { NoLoops = IGRAPH_NO_LOOPS, AllowLoops = IGRAPH_LOOPS };

enum EdgeOrder {
  EdgeById = IGRAPH_EDGEORDER_ID,
  EdgeBySourceId = IGRAPH_EDGEORDER_FROM,
  EdgeByTargetId = IGRAPH_EDGEORDER_TO
};

static inline int SetSeed(unsigned long int seed) {
  igraph_rng_t *rng = igraph_rng_default();
  int ret = igraph_rng_seed(rng, seed);
  return SafeCall(ret);
}

} // namespace igraph

#include "./vector.hpp"
#include "./vertex_selector.hpp"
#include "./vertex_iterator.hpp"
#include "./edge_selector.hpp"
#include "./edge_iterator.hpp"
#include "./graph.hpp"

#include "./vector_impl.hpp"
#include "./vertex_selector_impl.hpp"
#include "./vertex_iterator_impl.hpp"
#include "./edge_selector_impl.hpp"
#include "./edge_iterator_impl.hpp"
#include "./graph_impl.hpp"

#endif // IGRAPHPP_IGRAPH_HPP_
