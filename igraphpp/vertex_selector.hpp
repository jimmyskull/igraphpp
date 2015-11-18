#ifndef IGRAPHPP_VERTEX_SELECTOR_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <initializer_list>

#include <igraph.h>

#include "./util.hpp"

namespace igraph {

class Graph;

class VertexSelector {
public:
  ~VertexSelector();
  VertexSelector(const VertexSelector &vs);
  VertexSelector(std::initializer_list<double> list);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  VertexSelector(Iterator begin, Iterator end);
  VertexSelector(const VectorView &vector);

  bool is_all() const noexcept;
  int size(const Graph &graph) const noexcept;
  int type() const noexcept;

  /* Vertex selector constructors */
  static VertexSelector All();
  static VertexSelector Adjacent(int vid, NeighborMode mode = Out);
  static VertexSelector NonAdjacent(int vid, NeighborMode mode = Out);
  static VertexSelector None();
  static VertexSelector Single(int vid);
  // VectorView will not copy the contents of |vector|.
  static VertexSelector FromVector(const VectorView &vector);
  static VertexSelector Sequence(int from, int to);
  template <typename... Args, typename = typename std::enable_if<util::all_args(
                                  std::is_same<Args, int>::value...)>::type>
  static VertexSelector Small(Args... args);

  const igraph_vs_t &vs() const { return vs_; }
  igraph_vs_t *ptr() { return &vs_; }
  const igraph_vs_t *ptr() const { return &vs_; }

private:
  VertexSelector(const igraph_vs_t &vs);

  igraph_vs_t vs_;
};

template <typename Iterator, typename>
VertexSelector::VertexSelector(Iterator begin, Iterator end) {
  Vector vector(begin, end);
  SafeCall(igraph_vs_vector_copy(ptr(), vector.ptr()));
}

template <typename... Args, typename>
VertexSelector VertexSelector::Small(Args... args) {
  igraph_vs_t vs;
  SafeCall(igraph_vs_vector_small(&vs, args..., -1));
  return VertexSelector(vs);
}

} // namespace igraph

#endif // IGRAPHPP_VERTEX_SELECTOR_HPP_
