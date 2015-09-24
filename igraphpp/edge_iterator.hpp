#ifndef IGRAPHPP_EDGE_ITERATOR_HPP_
#define IGRAPHPP_EDGE_ITERATOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <iterator>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class EdgeSelector;

class EdgeIterator : public std::iterator<std::input_iterator_tag, long int> {
public:
  ~EdgeIterator();
  EdgeIterator(const Graph &graph, const EdgeSelector &es);
  EdgeIterator(const EdgeIterator &vit) = default;

  /* Stepping over the vertices */
  void next();
  bool at_end() const;
  long int size() const;
  void reset();

  /* Iterator */
  EdgeIterator &operator++();
  EdgeIterator operator++(int);
  bool operator==(const EdgeIterator &rhs) const;
  bool operator!=(const EdgeIterator &rhs) const;
  long int operator*();

  EdgeIterator begin() const;
  EdgeIterator end() const;

  const igraph_eit_t &eit() const { return eit_; }
  igraph_eit_t *ptr() { return &eit_; }
  const igraph_eit_t *ptr() const { return &eit_; }

private:
  EdgeIterator(const igraph_eit_t &eit);

  igraph_eit_t &eit() { return eit_; }

  igraph_eit_t eit_;
};

} // namespace igraph

#endif // IGRAPHPP_EDGE_ITERATOR_HPP_
