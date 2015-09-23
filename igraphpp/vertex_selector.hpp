#ifndef IGRAPHPP_VERTEX_SELECTOR_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <type_traits>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class VertexSelector {
public:
  ~VertexSelector() {
    if (view()) {
      igraph_vs_destroy(ptr());
    }
  }

  VertexSelector(const igraph_vs_t &vs, bool view) : view_(view) { vs_ = vs; }

  VertexSelector(const VertexSelector &vs) {
    SafeCall(igraph_vs_copy(ptr(), vs.ptr()));
  }

  bool is_all() const noexcept { return igraph_vs_is_all(ptr()); }

  // int size() const { return igraph_vs_size(ptr()); }

  int type() const noexcept { return igraph_vs_type(ptr()); }

  // static VertexSelector All() { return igraph_vss_all(); }

  // static VertexSelector None() { return igraph_vss_none(); }

  // static VertexSelector Adjacent(int vid, NeighborMode mode) {
  //   igraph_vs_t vs;
  //   SafeCall(igraph_vs_adj(&vs, vid, static_cast<igraph_neimode_t>(mode)));
  //   return VertexSelector(vs);
  // }

  // static VertexSelector NonAdjacent(int vid, NeighborMode mode) {
  //   igraph_vs_t vs;
  //   SafeCall(igraph_vs_nonadj(&vs, vid,
  //   static_cast<igraph_neimode_t>(mode)));
  //   return VertexSelector(vs);
  // }

  // static VertexSelector Single(int vid) { return igraph_vss_1(vid); }

  // static VertexSelector ViewVector(const VectorView &vector) {
  //   return igraph_vss_vector(vector.ptr());
  // }

  // static VertexSelector FromVector(const VectorView &vector) {
  //   igraph_vs_t vs;
  //   SafeCall(igraph_vs_vector_copy(&vs, vector.ptr()));
  //   return VertexSelector(vs);
  // }

  // constexpr bool all_args() { return true; }

  // template <typename... Tail> constexpr bool all_args(bool head, Tail...
  // tail) {
  //   return head && all_args(tail...);
  // }

  // template <typename... Args, typename>
  // static VertexSelector Small(const VectorView &vector, Args... args) {
  //   igraph_vs_t vs;
  //   SafeCall(igraph_vs_vector_small(&vs, args..., -1));
  //   return VertexSelector(vs);
  // }

  // static VertexSelector Sequence(int from, int to) { return
  // igraph_vss_seq(from, to); }

  igraph_vs_t *ptr() { return &vs_; }

  const igraph_vs_t *ptr() const { return &vs_; }

protected:
  igraph_vs_t vs_;

private:
  bool view() const { return view_; }

  bool view_ = false; /* Mark whether the instance is a view of a vector */
};

} // namespace igraph

#endif // IGRAPHPP_VERTEX_SELECTOR_HPP_
