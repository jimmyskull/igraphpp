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
  ~Graph() {
    if (VECTOR(graph_.from) != NULL)
      SafeCall(igraph_destroy(&graph_));
  }

  Graph(int vertices = 0, Directedness dir = Undirected) {
    SafeCall(igraph_empty(&graph_, vertices, dir));
  }

  Graph(const VectorView &vector, int vertices = 0, bool directed = false) {
    SafeCall(igraph_create(&graph_, vector.ptr(), vertices, directed));
  }

  Graph(const Graph &other) { SafeCall(igraph_copy(&graph_, &other.graph_)); }

  Graph(Graph &&other) {
    graph_ = other.graph_;
#ifndef NDEBUG
    std::memset(&other.graph_, 0, sizeof(other.graph_));
#endif
    VECTOR(other.graph_.from) = NULL;
  }

  Graph &operator=(const Graph &other) {
    if (this == &other)
      return *this;
    SafeCall(igraph_destroy(&graph_));
    SafeCall(igraph_copy(&graph_, &other.graph_));
    return *this;
  }

  Graph &operator=(Graph &&other) {
    graph_ = other.graph_;
#ifndef NDEBUG
    std::memset(&other.graph_, 0, sizeof(other.graph_));
#endif
    other.graph_.from.stor_begin = NULL;
    return *this;
  }

  int vcount() const noexcept { return igraph_vcount(&graph_); }

  int ecount() const noexcept { return igraph_ecount(&graph_); }

  bool is_directed() const noexcept { return igraph_is_directed(&graph_); }

  bool is_connected(Connectedness mode = WeaklyConnected) const {
    igraph_bool_t connected;
    igraph_connectedness_t m = static_cast<igraph_connectedness_t>(mode);
    int ret = igraph_is_connected(&graph_, &connected, m);
    SafeCall(ret);
    return static_cast<bool>(connected);
  }

  int diameter(Directedness directed = Directed,
               bool unconnected = true) const {
    int length = 0;
    int ret = igraph_diameter(&graph_, &length, NULL, NULL, NULL, directed,
                              unconnected);
    SafeCall(ret);
    return length;
  }

  Vector degrees(NeighborMode mode = All, Loops loops = AllowLoops) const {
    Vector v;
    bool lps = (loops == AllowLoops);
    SafeCall(igraph_degree(&graph_, v.ptr(), igraph_vss_all(),
                           static_cast<igraph_neimode_t>(mode), lps));
    return std::move(v);
  }

  Graph &AddEdge(int from, int to) {
    SafeCall(igraph_add_edge(&graph_, from, to));
    return *this;
  }

  Graph &AddEdges(const Vector &edges) {
    SafeCall(igraph_add_edges(&graph_, edges.ptr(), NULL));
    return *this;
  }

  double AveragePathLength(Directedness directed = Directed,
                           bool unconnected = true) const {
    double ret;
    SafeCall(igraph_average_path_length(&graph_, &ret, directed, unconnected));
    return ret;
  }

  static Graph Lattice(const Vector &dimension, int nei = 1,
                       Directedness dir = Undirected,
                       Mutuality mutual = NotMutual,
                       Periodicity periodicity = NotPeriodic) {
    igraph_t graph;
    int ret =
        igraph_lattice(&graph, dimension.ptr(), nei, dir, mutual, periodicity);
    SafeCall(ret);
    return Graph(graph);
  }

  static Graph ErdosRenyiGame(int vertices, double prob,
                              Directedness dir = Undirected,
                              Loops loops = NoLoops) {
    igraph_t graph;
    int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, vertices,
                                      prob, dir, loops);
    SafeCall(ret);
    return Graph(graph);
  }

  static Graph ErdosRenyiGame(int vertices, int edges,
                              Directedness dir = Undirected,
                              Loops loops = NoLoops) {
    igraph_t graph;
    int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, vertices,
                                      edges, dir, loops);
    SafeCall(ret);
    return Graph(graph);
  }

protected:
  Graph(const igraph_t &graph) { SafeCall(igraph_copy(&graph_, &graph)); }

private:
  igraph_t graph_;
};

} // namespace igraph

#endif // IGRAPHPP_GRAPH_HPP_
