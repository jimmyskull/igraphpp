#ifndef IGRAPHPP_VERTEX_ITERATOR_HPP_
#define IGRAPHPP_VERTEX_ITERATOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <iterator>

#include <igraph/igraph.h>

#include "./exception.hpp"

namespace igraph {

class VertexSelector;

class VertexIterator : public std::iterator<std::input_iterator_tag, long int> {
 public:
  ~VertexIterator();
  VertexIterator(const Graph &graph, const VertexSelector &vs);
  VertexIterator(const VertexIterator &vit) = default;

  /* Stepping over the vertices */
  void next();
  bool at_end() const;
  long int size() const;
  void reset();

  /* Iterator */
  VertexIterator &operator++();
  VertexIterator operator++(int);
  bool operator==(const VertexIterator &rhs) const;
  bool operator!=(const VertexIterator &rhs) const;
  long int operator*();

  VertexIterator begin() const;
  VertexIterator end() const;

  const igraph_vit_t &vit() const { return vit_; }
  igraph_vit_t *ptr() { return &vit_; }
  const igraph_vit_t *ptr() const { return &vit_; }

 private:
  VertexIterator(const igraph_vit_t &vit);

  igraph_vit_t &vit() { return vit_; }

  igraph_vit_t vit_;
};

}  // namespace igraph

#endif  // IGRAPHPP_VERTEX_ITERATOR_HPP_
