#ifndef IGRAPHPP_EDGE_ITERATOR_IMPL_HPP_
#define IGRAPHPP_EDGE_ITERATOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline EdgeIterator::~EdgeIterator() { igraph_eit_destroy(ptr()); }
inline EdgeIterator::EdgeIterator(const Graph &graph, const EdgeSelector &es) {
  SafeCall(igraph_eit_create(graph.ptr(), es.es(), ptr()));
}

inline void EdgeIterator::next() { IGRAPH_EIT_NEXT(eit()); }
inline bool EdgeIterator::at_end() const { return IGRAPH_EIT_END(eit()); }
inline long int EdgeIterator::size() const { return IGRAPH_EIT_SIZE(eit()); }
inline void EdgeIterator::reset() { IGRAPH_EIT_RESET(eit()); }

inline EdgeIterator &EdgeIterator::operator++() {
  next();
  return *this;
}
inline EdgeIterator EdgeIterator::operator++(int) {
  EdgeIterator result(*this);
  next();
  return result;
}
inline bool EdgeIterator::operator==(const EdgeIterator &rhs) const {
  return eit().pos == rhs.eit().pos;
}
inline bool EdgeIterator::operator!=(const EdgeIterator &rhs) const {
  return !(*this == rhs);
}
inline long int EdgeIterator::operator*() { return IGRAPH_EIT_GET(eit()); }

inline EdgeIterator EdgeIterator::begin() const {
  igraph_eit_t copy = eit_;
  IGRAPH_VIT_RESET(copy);
  return EdgeIterator(copy);
}
inline EdgeIterator EdgeIterator::end() const {
  igraph_eit_t copy = eit_;
  copy.pos = copy.end;
  return EdgeIterator(copy);
}

inline EdgeIterator::EdgeIterator(const igraph_eit_t &eit) : eit_(eit) {}

} // namespace igraph

#endif // IGRAPHPP_EDGE_ITERATOR_IMPL_HPP_
