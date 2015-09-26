#ifndef IGRAPHPP_GRAPH_HPP_
#define IGRAPHPP_GRAPH_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <initializer_list>

#include <igraph.h>

#include "./exception.hpp"
#include "./util.hpp"

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
  Graph(std::initializer_list<double> edges, long int vertices = 0,
        Directedness dir = Undirected);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  Graph(Iterator edges_begin, Iterator edges_end, long int vertices = 0,
        Directedness dir = Undirected);
  Graph(const Graph &other);
  Graph(Graph &&other);

  Graph &operator=(const Graph &other);
  Graph &operator=(Graph &&other);

  /* Basic query operations */
  int vcount() const noexcept;
  int ecount() const noexcept;
  Edge edge(int eid) const;
  int eid(int from, int to, Directedness dir = Undirected) const;
  // Edge eids between pairs of vertices
  Vector pairs_eids(const VectorView &pairs, Directedness dir = Directed,
                    MultiEdges multiedges = IgnoreMultiEdges) const;
  Vector pairs_eids(std::initializer_list<double> pairs,
                    Directedness dir = Directed,
                    MultiEdges multiedges = IgnoreMultiEdges) const;
  // Edges composing a path
  Vector path_eids(const VectorView &path, Directedness dir = Directed,
                   MultiEdges multiedges = IgnoreMultiEdges) const;
  Vector path_eids(std::initializer_list<double> path,
                   Directedness dir = Directed,
                   MultiEdges multiedges = IgnoreMultiEdges) const;
  // Concatenation of |pairs_eids()| and |path_eids()|
  Vector eids(const VectorView &pairs, const VectorView &path,
              Directedness dir = Directed,
              MultiEdges multiedges = IgnoreMultiEdges) const;
  Vector eids(std::initializer_list<double> pairs,
              std::initializer_list<double> path, Directedness dir = Directed,
              MultiEdges multiedges = IgnoreMultiEdges) const;
  Vector neighbors(int vertex, NeighborMode mode = Out) const;
  Vector incident(int vertex, NeighborMode mode = Out) const;
  bool is_directed() const noexcept;
  Vector degree(const VertexSelector &vids, NeighborMode mode = Out,
                Loops loops = NoLoops) const;
  Vector degree(std::initializer_list<double> vids, NeighborMode mode = Out,
                Loops loops = NoLoops) const;

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

  igraph_t *ptr() { return &graph_; }
  const igraph_t *ptr() const { return &graph_; }

  Graph(const igraph_t *graph) {
    graph_ = *graph;
    disown();
  }

protected:
  Graph(const igraph_t &graph);

private:
  void disown() { owner_ = false; }
  bool owner() const { return owner_; }

  igraph_t graph_;
  bool owner_ = true;
};

} // namespace igraph

#endif // IGRAPHPP_GRAPH_HPP_
