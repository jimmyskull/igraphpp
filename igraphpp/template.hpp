#ifndef IGRAPHPP_VERTEX_ITERATOR_HPP_
#define IGRAPHPP_VERTEX_ITERATOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

class VertexIterator {
 public:
  ~VertexIterator();

  VertexIterator(const VertexIterator &vs);

  igraph_vit_t *ptr() { return &vit_; }

  const igraph_vit_t *ptr() const { return &vit_; }

 private:
  // VertexIterator(const igraph_vit_t &vs);

  igraph_vit_t vit_;
};

}  // namespace igraph

#endif  // IGRAPHPP_VERTEX_ITERATOR_HPP_
