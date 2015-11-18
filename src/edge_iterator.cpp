#include "../igraphpp/igraph.hpp"

#include <igraph.h>

namespace igraph {

EdgeIterator::~EdgeIterator() { igraph_eit_destroy(ptr()); }
EdgeIterator::EdgeIterator(const Graph &graph, const EdgeSelector &es) {
  SafeCall(igraph_eit_create(graph.ptr(), es.es(), ptr()));
}

void EdgeIterator::next() { IGRAPH_EIT_NEXT(eit()); }
bool EdgeIterator::at_end() const { return IGRAPH_EIT_END(eit()); }
long int EdgeIterator::size() const { return IGRAPH_EIT_SIZE(eit()); }
void EdgeIterator::reset() { IGRAPH_EIT_RESET(eit()); }

EdgeIterator &EdgeIterator::operator++() {
  next();
  return *this;
}
EdgeIterator EdgeIterator::operator++(int) {
  EdgeIterator result(*this);
  next();
  return result;
}
bool EdgeIterator::operator==(const EdgeIterator &rhs) const {
  return eit().pos == rhs.eit().pos;
}
bool EdgeIterator::operator!=(const EdgeIterator &rhs) const {
  return !(*this == rhs);
}
long int EdgeIterator::operator*() { return IGRAPH_EIT_GET(eit()); }

EdgeIterator EdgeIterator::begin() const {
  igraph_eit_t copy = eit_;
  IGRAPH_VIT_RESET(copy);
  return EdgeIterator(copy);
}
EdgeIterator EdgeIterator::end() const {
  igraph_eit_t copy = eit_;
  copy.pos = copy.end;
  return EdgeIterator(copy);
}

EdgeIterator::EdgeIterator(const igraph_eit_t &eit) : eit_(eit) {}

}  // namespace igraph
