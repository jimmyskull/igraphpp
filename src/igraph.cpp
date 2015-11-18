#include "../igraphpp/igraph.hpp"

#include "../igraphpp/exception.hpp"

namespace igraph {

int SetSeed(unsigned long int seed) {
  igraph_rng_t *rng = igraph_rng_default();
  int ret = igraph_rng_seed(rng, seed);
  return SafeCall(ret);
}

}  // namespace igraph
