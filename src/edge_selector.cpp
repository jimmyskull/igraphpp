#include "../igraphpp/igraph.hpp"

#include <type_traits>

#include <igraph.h>

namespace igraph {

EdgeSelector::~EdgeSelector() { igraph_es_destroy(ptr()); }

EdgeSelector::EdgeSelector(const EdgeSelector &es) {
  SafeCall(igraph_es_copy(ptr(), es.ptr()));
}
EdgeSelector::EdgeSelector(std::initializer_list<double> list) {
  Vector vector(list);
  SafeCall(igraph_es_vector_copy(ptr(), vector.ptr()));
}
EdgeSelector::EdgeSelector(const VectorView &vector) {
  SafeCall(igraph_es_vector_copy(ptr(), vector.ptr()));
}

bool EdgeSelector::is_all() const noexcept { return igraph_es_is_all(ptr()); }

int EdgeSelector::size(const Graph &graph) const noexcept {
  int result;
  SafeCall(igraph_es_size(graph.ptr(), ptr(), &result));
  return result;
}

int EdgeSelector::type() const noexcept { return igraph_es_type(ptr()); }

/* Vertex selector constructors */
EdgeSelector EdgeSelector::All(EdgeOrder order) {
  static EdgeSelector instance(
      igraph_ess_all(static_cast<igraph_edgeorder_type_t>(order)));
  return instance;
}
EdgeSelector EdgeSelector::Incident(int vid, NeighborMode mode) {
  igraph_es_t es;
  SafeCall(igraph_es_incident(&es, vid, static_cast<igraph_neimode_t>(mode)));
  return EdgeSelector(es);
}
EdgeSelector EdgeSelector::None() {
  static EdgeSelector instance(igraph_ess_none());
  return instance;
}
EdgeSelector EdgeSelector::Single(int eid) {
  static EdgeSelector instance(igraph_ess_1(eid));
  return instance;
}
EdgeSelector EdgeSelector::FromVector(const VectorView &vector) {
  static EdgeSelector instance(igraph_ess_vector(vector.ptr()));
  return instance;
}
EdgeSelector EdgeSelector::Sequence(int from, int to) {
  static EdgeSelector instance(igraph_ess_seq(from, to));
  return instance;
}
EdgeSelector EdgeSelector::Pairs(const Vector &edge_vector, Directedness dir) {
  igraph_es_t es;
  SafeCall(igraph_es_pairs(&es, edge_vector.ptr(), dir));
  return EdgeSelector(es);
}

EdgeSelector::EdgeSelector(const igraph_es_t &es) { es_ = es; }

}  // namespace igraph
