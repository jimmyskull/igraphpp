#ifndef IGRAPHPP_GRAPH_HPP_
#define IGRAPHPP_GRAPH_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class Graph {
public:
  ~Graph();
  Graph(int vertices = 0, Directedness dir = Undirected);

  Graph(const VectorView &vector, int vertices = 0, bool directed = false);

  Graph(const Graph &other);

  Graph(Graph &&other);

  Graph &operator=(const Graph &other);

  Graph &operator=(Graph &&other);

  int vcount() const noexcept;

  int ecount() const noexcept;

  bool is_directed() const noexcept;

  bool is_connected(Connectedness mode = WeaklyConnected) const;

  int diameter(Directedness directed = Directed, bool unconnected = true) const;

  Vector degrees(NeighborMode mode = All, Loops loops = AllowLoops) const;

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

  void set_graph_destroy(bool value) { graph_destroy_ = value; }

  igraph_t *ptr() { return &graph_; }

  const igraph_t *ptr() const { return &graph_; }

  Graph(const igraph_t *graph) {
    graph_ = *graph;
    set_graph_destroy(false);
  }

protected:
  Graph(const igraph_t &graph);

private:
  igraph_t graph_;
  bool graph_destroy_ = true;
};

} // namespace igraph

#endif // IGRAPHPP_GRAPH_HPP_
