#include "../igraphpp/igraph.hpp"

#include <type_traits>

#include <igraph.h>

namespace igraph {

VertexSelector::~VertexSelector() { igraph_vs_destroy(ptr()); }
VertexSelector::VertexSelector(const VertexSelector &vs) {
  SafeCall(igraph_vs_copy(ptr(), vs.ptr()));
}
VertexSelector::VertexSelector(std::initializer_list<double> list) {
  Vector vector(list);
  SafeCall(igraph_vs_vector_copy(ptr(), vector.ptr()));
}
VertexSelector::VertexSelector(const VectorView &vector) {
  SafeCall(igraph_vs_vector_copy(ptr(), vector.ptr()));
}

bool VertexSelector::is_all() const noexcept { return igraph_vs_is_all(ptr()); }
int VertexSelector::size(const Graph &graph) const noexcept {
  int result;
  SafeCall(igraph_vs_size(graph.ptr(), ptr(), &result));
  return result;
}
int VertexSelector::type() const noexcept { return igraph_vs_type(ptr()); }

/* Vertex selector constructors */
VertexSelector VertexSelector::All() {
  static VertexSelector instance(igraph_vss_all());
  return instance;
}
VertexSelector VertexSelector::Adjacent(int vid, Mode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_adj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}
VertexSelector VertexSelector::NonAdjacent(int vid, Mode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_nonadj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}
VertexSelector VertexSelector::None() {
  static VertexSelector instance(igraph_vss_none());
  return instance;
}
VertexSelector VertexSelector::Single(int vid) {
  return VertexSelector(igraph_vss_1(vid));
}
VertexSelector VertexSelector::FromVector(const VectorView &vector) {
  return VertexSelector(igraph_vss_vector(vector.ptr()));
}
VertexSelector VertexSelector::Sequence(int from, int to) {
  return VertexSelector(igraph_vss_seq(from, to));
}

VertexSelector::VertexSelector(const igraph_vs_t &vs) { vs_ = vs; }

}  // namespace igraph
