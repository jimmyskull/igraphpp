#ifndef IGRAPHPP_VERTEX_SELECTOR_HPP_
#define IGRAPHPP_VERTEX_SELECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <type_traits>

#include <igraph.h>

namespace igraph {

namespace {

constexpr bool all_args() { return true; }

template <typename... Tail> constexpr bool all_args(bool head, Tail... tail) {
  return head && all_args(tail...);
}

} // namespace

class Graph;

class VertexSelector {
public:
  ~VertexSelector();
  VertexSelector(const VertexSelector &vs);

  bool is_all() const noexcept;
  int size(const Graph &graph) const noexcept;
  int type() const noexcept;

  /* Vertex selector constructors */
  static VertexSelector All();
  static VertexSelector Adjacent(int vid, NeighborMode mode = Out);
  static VertexSelector NonAdjacent(int vid, NeighborMode mode = Out);
  static VertexSelector None();
  static VertexSelector Single(int vid);
  static VertexSelector FromVector(const VectorView &vector);
  static VertexSelector Sequence(int from, int to);

  template <typename... Args, typename = std::enable_if_t<
                                  all_args(std::is_same<Args, int>::value...)>>
  static VertexSelector Small(Args... args);

  const igraph_vs_t &vs() const { return vs_; }
  igraph_vs_t *ptr() { return &vs_; }
  const igraph_vs_t *ptr() const { return &vs_; }

private:
  VertexSelector(const igraph_vs_t &vs);

  igraph_vs_t vs_;
};

} // namespace igraph

#endif // IGRAPHPP_VERTEX_SELECTOR_HPP_
