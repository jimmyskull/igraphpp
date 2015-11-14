#ifndef IGRAPHPP_EDGE_SELECTOR_HPP_
#define IGRAPHPP_EDGE_SELECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <type_traits>

#include <igraph.h>

#include "./util.hpp"

namespace igraph {

class Graph;

class Vector;

class EdgeSelector {
 public:
  ~EdgeSelector();
  EdgeSelector(const EdgeSelector &es);
  EdgeSelector(std::initializer_list<double> list);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  EdgeSelector(Iterator begin, Iterator end);
  EdgeSelector(const VectorView &vector);

  bool is_all() const noexcept;
  int size(const Graph &graph) const noexcept;
  int type() const noexcept;

  /* Vertex selector constructors */
  static EdgeSelector All(EdgeOrder order = EdgeById);
  static EdgeSelector Incident(int vid, NeighborMode mode = Out);
  static EdgeSelector None();
  static EdgeSelector Single(int eid);
  static EdgeSelector FromVector(const VectorView &vector);
  static EdgeSelector Sequence(int from, int to);
  static EdgeSelector Pairs(const Vector &edge_vector,
                            Directedness dir = Directed);
  template <typename... Args, typename = typename std::enable_if<util::all_args(
                                  std::is_same<Args, int>::value...)>::type>
  static EdgeSelector Pairs(Directedness dir, Args... args);

  const igraph_es_t &es() const { return es_; }
  igraph_es_t *ptr() { return &es_; }
  const igraph_es_t *ptr() const { return &es_; }

 protected:
  EdgeSelector(const igraph_es_t &es);

 private:
  igraph_es_t es_;
};

template <typename Iterator, typename>
EdgeSelector::EdgeSelector(Iterator begin, Iterator end) {
  Vector vector(begin, end);
  SafeCall(igraph_es_vector_copy(ptr(), vector.ptr()));
}

template <typename... Args, typename>
EdgeSelector EdgeSelector::Pairs(Directedness dir, Args... args) {
  igraph_es_t es;
  SafeCall(igraph_es_pairs_small(&es, dir, args..., -1));
  return EdgeSelector(es);
}

}  // namespace igraph

#endif  // IGRAPHPP_EDGE_SELECTOR_HPP_
