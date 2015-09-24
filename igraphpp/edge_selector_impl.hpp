#ifndef IGRAPHPP_EDGE_SELECTOR_IMPL_HPP_
#define IGRAPHPP_EDGE_SELECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <type_traits>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline EdgeSelector::~EdgeSelector() { igraph_es_destroy(ptr()); }

inline EdgeSelector::EdgeSelector(const EdgeSelector &es) {
  SafeCall(igraph_es_copy(ptr(), es.ptr()));
}

inline bool EdgeSelector::is_all() const noexcept {
  return igraph_es_is_all(ptr());
}

inline int EdgeSelector::size(const Graph &graph) const noexcept {
  int result;
  SafeCall(igraph_es_size(graph.ptr(), ptr(), &result));
  return result;
}

inline int EdgeSelector::type() const noexcept { return igraph_es_type(ptr()); }

/* Vertex selector constructors */
inline EdgeSelector EdgeSelector::All(EdgeOrder order) {
  static EdgeSelector instance(
      igraph_ess_all(static_cast<igraph_edgeorder_type_t>(order)));
  return instance;
}
inline EdgeSelector EdgeSelector::Incident(int vid, NeighborMode mode) {
  igraph_es_t es;
  SafeCall(igraph_es_incident(&es, vid, static_cast<igraph_neimode_t>(mode)));
  return EdgeSelector(es);
}
inline EdgeSelector EdgeSelector::None() {
  static EdgeSelector instance(igraph_ess_none());
  return instance;
}
inline EdgeSelector EdgeSelector::Single(int eid) {
  static EdgeSelector instance(igraph_ess_1(eid));
  return instance;
}
inline EdgeSelector EdgeSelector::FromVector(const VectorView &vector) {
  static EdgeSelector instance(igraph_ess_vector(vector.ptr()));
  return instance;
}
inline EdgeSelector EdgeSelector::Sequence(int from, int to) {
  static EdgeSelector instance(igraph_ess_seq(from, to));
  return instance;
}
inline EdgeSelector EdgeSelector::Pairs(const Vector &edge_vector,
                                        Directedness dir) {
  igraph_es_t es;
  SafeCall(igraph_es_pairs(&es, edge_vector.ptr(), dir));
  return EdgeSelector(es);
}
template <typename... Args, typename>
inline EdgeSelector EdgeSelector::Pairs(Directedness dir, Args... args) {
  igraph_es_t es;
  SafeCall(igraph_es_pairs_small(&es, dir, args..., -1));
  return EdgeSelector(es);
}

inline EdgeSelector::EdgeSelector(const igraph_es_t &es) { es_ = es; }

} // namespace igraph

#endif // IGRAPHPP_EDGE_SELECTOR_IMPL_HPP_
