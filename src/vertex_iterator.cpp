#include "../igraphpp/igraph.hpp"

#include <igraph.h>

namespace igraph {

VertexIterator::~VertexIterator() { igraph_vit_destroy(ptr()); }
VertexIterator::VertexIterator(const Graph &graph, const VertexSelector &vs) {
  SafeCall(igraph_vit_create(graph.ptr(), vs.vs(), ptr()));
}

void VertexIterator::next() { IGRAPH_VIT_NEXT(vit()); }
bool VertexIterator::at_end() const { return IGRAPH_VIT_END(vit()); }
long int VertexIterator::size() const { return IGRAPH_VIT_SIZE(vit()); }
void VertexIterator::reset() { IGRAPH_VIT_RESET(vit()); }

VertexIterator &VertexIterator::operator++() {
  next();
  return *this;
}
VertexIterator VertexIterator::operator++(int) {
  VertexIterator result(*this);
  next();
  return result;
}
bool VertexIterator::operator==(const VertexIterator &rhs) const {
  return vit().pos == rhs.vit().pos;
}
bool VertexIterator::operator!=(const VertexIterator &rhs) const {
  return !(*this == rhs);
}
long int VertexIterator::operator*() { return IGRAPH_VIT_GET(vit()); }

VertexIterator VertexIterator::begin() const {
  igraph_vit_t copy = vit_;
  IGRAPH_VIT_RESET(copy);
  return VertexIterator(copy);
}
VertexIterator VertexIterator::end() const {
  igraph_vit_t copy = vit_;
  copy.pos = copy.end;
  return VertexIterator(copy);
}

VertexIterator::VertexIterator(const igraph_vit_t &vit) : vit_(vit) {}

}  // namespace igraph
