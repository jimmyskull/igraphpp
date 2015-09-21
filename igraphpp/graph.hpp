#ifndef IGRAPHPP_GRAPH_HPP_
#define IGRAPHPP_GRAPH_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class Graph {
public:
  Graph(int vertices = 0, Directedness dir = Undirected);

  Graph(const VectorView &vector, int vertices = 0, bool directed = false);

  ~Graph();

  Graph(const Graph &other);

  Graph(Graph &&other);

  Graph &operator=(const Graph &other);

  Graph &operator=(Graph &&other);

  int vcount() const noexcept;

  int ecount() const noexcept;

  bool is_directed() const noexcept;

  bool is_connected(Connectedness mode = WeaklyConnected) const;

  int diameter(Directedness directed = Directed, bool unconnected = true) const;

  Vector degrees(DegreeMode mode = All, Loops loops = AllowLoops) const;

  Graph &AddEdge(int from, int to);

  Graph &AddEdges(const Vector &edges);

  double AveragePathLength(Directedness directed = Directed,
                           bool unconnected = true) const;

  static Graph Lattice(const Vector &dimension, int nei = 1,
                       Directedness dir = Undirected,
                       Mutuality mutual = NotMutual,
                       Periodicity periodicity = NotPeriodic);

  static Graph ErdosRenyiGame(int vertices, double prob,
                              Directedness dir = Undirected,
                              Loops loops = NoLoops);

  static Graph ErdosRenyiGame(int vertices, int edges,
                              Directedness dir = Undirected,
                              Loops loops = NoLoops);

protected:
  Graph(const igraph_t &graph);

private:
  igraph_t graph_;
};

} // namespace igraph

#include "./graph_impl.hpp"

#endif // IGRAPHPP_GRAPH_HPP_
