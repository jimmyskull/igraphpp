#ifndef IGRAPHPP_VERTEX_SELECTOR_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class VertexSelector {
public:
  ~VertexSelector();
  VertexSelector(const igraph_vs_t &vs, bool view = false);
  VertexSelector(const VertexSelector &vs);

  /* Generic vertex selector operations */
  bool is_all() const noexcept;
  // int size() const;
  int type() const noexcept;

  /* Vertex selector constructors */
  static VertexSelector All();
  static VertexSelector None();
  static VertexSelector Adjacent(int vid, NeighborMode mode = Out);
  static VertexSelector NonAdjacent(int vid, NeighborMode mode = Out);
  static VertexSelector Single(int vid);
  /* A view of |vector|, sharing its memory. Will not release memory. */
  static VertexSelector ViewVector(const VectorView &vector);
  /* Makes a copy of |vector|. */
  static VertexSelector FromVector(const VectorView &vector);
  template <typename... Args>
  static VertexSelector Small(const VectorView &vector, Args &&... args);
  static VertexSelector Sequence(int from, int to);

  igraph_vs_t *ptr();
  const igraph_vs_t *ptr() const;

protected:
  igraph_vs_t vs_;

private:
  bool view() const { return view_; }

  bool view_ = false; /* Mark whether the instance is a view of a vector */
};

} // namespace igraph

#include "./vertex_selector_impl.hpp"

#endif // IGRAPHPP_VERTEX_SELECTOR_HPP_
