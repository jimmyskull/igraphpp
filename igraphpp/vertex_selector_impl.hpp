#ifndef IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

VertexSelector::~VertexSelector() {
  if (view()) {
    igraph_vs_destroy(ptr());
  }
}

VertexSelector::VertexSelector(const igraph_vs_t &vs, bool view) : view_(view) {
  vs_ = vs;
}

VertexSelector::VertexSelector(const VertexSelector &vs) {
  SafeCall(igraph_vs_copy(ptr(), vs.ptr()));
}

bool VertexSelector::is_all() const noexcept { return igraph_vs_is_all(ptr()); }

// int VertexSelector::size() const { return igraph_vs_size(ptr()); }

int VertexSelector::type() const noexcept { return igraph_vs_type(ptr()); }

VertexSelector VertexSelector::All() { return igraph_vss_all(); }

VertexSelector VertexSelector::None() { return igraph_vss_none(); }

VertexSelector VertexSelector::Adjacent(int vid, NeighborMode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_adj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}

VertexSelector VertexSelector::NonAdjacent(int vid, NeighborMode mode) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_nonadj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  return VertexSelector(vs);
}

VertexSelector VertexSelector::Single(int vid) { return igraph_vss_1(vid); }

VertexSelector VertexSelector::ViewVector(const VectorView &vector) {
  return igraph_vss_vector(vector.ptr());
}

VertexSelector VertexSelector::FromVector(const VectorView &vector) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_vector_copy(&vs, vector.ptr()));
  return VertexSelector(vs);
}

template <typename... Args>
VertexSelector VertexSelector::Small(const VectorView &vector,
                                     Args &&... args) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_vector_small(&vs, args..., -1));
  return VertexSelector(vs);
}

VertexSelector VertexSelector::Sequence(int from, int to) {
  return igraph_vss_seq(from, to);
}

igraph_vs_t *VertexSelector::ptr() { return &vs_; }

const igraph_vs_t *VertexSelector::ptr() const { return &vs_; }

} // namespace igraph

#endif // IGRAPHPP_VERTEX_SELECTOR_IMPL_HPP_
