#ifndef IGRAPHPP_VERTEX_ITERATOR_IMPL_HPP_
#define IGRAPHPP_VERTEX_ITERATOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline VertexIterator::~VertexIterator() { igraph_vit_destroy(ptr()); }
inline VertexIterator::VertexIterator(const Graph &graph,
                                      const VertexSelector &vs) {
  SafeCall(igraph_vit_create(graph.ptr(), vs.vs(), ptr()));
}

inline void VertexIterator::next() { IGRAPH_VIT_NEXT(vit()); }
inline bool VertexIterator::at_end() const { return IGRAPH_VIT_END(vit()); }
inline long int VertexIterator::size() const { return IGRAPH_VIT_SIZE(vit()); }
inline void VertexIterator::reset() { IGRAPH_VIT_RESET(vit()); }

inline VertexIterator &VertexIterator::operator++() {
  next();
  return *this;
}
inline VertexIterator VertexIterator::operator++(int) {
  VertexIterator result(*this);
  next();
  return result;
}
inline bool VertexIterator::operator==(const VertexIterator &rhs) const {
  return vit().pos == rhs.vit().pos;
}
inline bool VertexIterator::operator!=(const VertexIterator &rhs) const {
  return !(*this == rhs);
}
inline long int VertexIterator::operator*() { return IGRAPH_VIT_GET(vit()); }

inline VertexIterator VertexIterator::begin() const {
  igraph_vit_t copy = vit_;
  IGRAPH_VIT_RESET(copy);
  return VertexIterator(copy);
}
inline VertexIterator VertexIterator::end() const {
  igraph_vit_t copy = vit_;
  copy.pos = copy.end;
  return VertexIterator(copy);
}

inline VertexIterator::VertexIterator(const igraph_vit_t &vit) : vit_(vit) {}

} // namespace igraph

#endif // IGRAPHPP_VERTEX_ITERATOR_IMPL_HPP_
