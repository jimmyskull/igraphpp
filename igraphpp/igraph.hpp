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

enum Loops { NoLoops = IGRAPH_NO_LOOPS, AllowLoops = IGRAPH_LOOPS };

static inline int SetSeed(unsigned long int seed) {
  igraph_rng_t *rng = igraph_rng_default();
  int ret = igraph_rng_seed(rng, seed);
  return SafeCall(ret);
}

} // namespace igraph

#include "./vector.hpp"
#include "./graph.hpp"

#endif // IGRAPHPP_IGRAPH_HPP_
