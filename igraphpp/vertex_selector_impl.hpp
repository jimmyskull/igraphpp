#ifndef IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <type_traits>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline VertexSelector::~VertexSelector() { igraph_vs_destroy(ptr()); }

inline VertexSelector::VertexSelector(const VertexSelector &vs) {
  SafeCall(igraph_vs_copy(ptr(), vs.ptr()));
}

inline bool VertexSelector::is_all() const noexcept {
  return igraph_vs_is_all(ptr());
}

inline int VertexSelector::size(const Graph &graph) const noexcept {
  int result;
  SafeCall(igraph_vs_size(graph.ptr(), ptr(), &result));
  return result;
}

inline int VertexSelector::type() const noexcept {
  return igraph_vs_type(ptr());
}

/* Vertex selector constructors */
inline VertexSelector VertexSelector::All() {
  static VertexSelector instance(igraph_vss_all());
  return instance;
}
inline VertexSelector VertexSelector::Adjacent(int vid, NeighborMode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_adj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}
inline VertexSelector VertexSelector::NonAdjacent(int vid, NeighborMode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_nonadj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}
inline VertexSelector VertexSelector::None() {
  static VertexSelector instance(igraph_vss_none());
  return instance;
}
inline VertexSelector VertexSelector::Single(int vid) {
  static VertexSelector instance(igraph_vss_1(vid));
  return instance;
}
inline VertexSelector VertexSelector::FromVector(const VectorView &vector) {
  static VertexSelector instance(igraph_vss_vector(vector.ptr()));
  return instance;
}
inline VertexSelector VertexSelector::Sequence(int from, int to) {
  static VertexSelector instance(igraph_vss_seq(from, to));
  return instance;
}
template <typename... Args, typename>
inline VertexSelector VertexSelector::Small(Args... args) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_vector_small(&vs, args..., -1));
  return VertexSelector(vs);
}

inline VertexSelector::VertexSelector(const igraph_vs_t &vs) { vs_ = vs; }

} // namespace igraph

#endif // IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_
