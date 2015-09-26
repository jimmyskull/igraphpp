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

class Vector;
class VertexSelector;

class Graph {
public:
  /* Constructors and Destructors */
  ~Graph();
  Graph(long int vertices = 0, Directedness dir = Undirected);
  // TODO: Empty graph with attributes
  Graph(const VectorView &edges, long int vertices = 0,
        Directedness dir = Undirected);
  Graph(const Graph &other);
  Graph(Graph &&other);

  Graph &operator=(const Graph &other);
  Graph &operator=(Graph &&other);

  /* Basic query operations */
  int vcount() const noexcept;
  int ecount() const noexcept;
  // Edge edge(int eid) const noexcept;
  // int eid(int from, int to, Directedness directed = Undirected) const;
  // Vector pairs_eids(const Vector &pairs, Directedness directed = Undirected,
  //                   bool fail_if_not_connected = false,
  //                   bool multiedges = false) const;
  // Vector path_eids(const Vector &path, Directedness directed = Undirected,
  //                  bool fail_if_not_connected = false,
  //                  bool multiedges = false) const;
  // Vector eids(const Vector &pairs, const Vector &path,
  //             Directedness directed = Undirected,
  //             bool fail_if_not_connected = false,
  //             bool multiedges = false) const;
  // Vector neighbors(int vertex, Directedness directed = Undirected) const;
  // Vector incident(int vertex, Directedness directed = Undirected) const;
  bool is_directed() const noexcept;
  // Vector degree(const VertexSelector &vids,
  //               NeighborMode mode = All, Loops loops = AllowLoops) const;

  /* Adding and deleting vertices and edges */
  Graph &AddEdge(int from, int to);
  Graph &AddEdges(const Vector &edges);

  bool is_connected(Connectedness mode = WeaklyConnected) const;

  int diameter(Directedness directed = Directed, bool unconnected = true) const;

  double AveragePathLength(Directedness directed = Directed,
                           bool unconnected = true) const;

  static Graph Empty(int vertices = 0, Directedness dir = Undirected);
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
  void disown() { VECTOR(graph_.from) = NULL; }
  bool owner() const { return graph_destroy_ && VECTOR(graph_.from) != NULL; }

  igraph_t graph_;
  bool graph_destroy_ = true;
};

} // namespace igraph

#endif // IGRAPHPP_GRAPH_HPP_
