#ifndef IGRAPHPP_GRAPH_IMPL_HPP_
#define IGRAPHPP_GRAPH_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline Graph::~Graph() {
  if (graph_destroy_ && VECTOR(graph_.from) != NULL)
    SafeCall(igraph_destroy(&graph_));
}

inline Graph::Graph(int vertices, Directedness dir) {
  SafeCall(igraph_empty(&graph_, vertices, dir));
}

inline Graph::Graph(const VectorView &vector, int vertices, bool directed) {
  SafeCall(igraph_create(&graph_, vector.ptr(), vertices, directed));
}

inline Graph::Graph(const Graph &other) {
  SafeCall(igraph_copy(&graph_, &other.graph_));
}

inline Graph::Graph(Graph &&other) {
  graph_ = other.graph_;
#ifndef NDEBUG
  std::memset(&other.graph_, 0, sizeof(other.graph_));
#endif
  VECTOR(other.graph_.from) = NULL;
}

inline Graph &Graph::operator=(const Graph &other) {
  if (this == &other)
    return *this;
  SafeCall(igraph_destroy(&graph_));
  SafeCall(igraph_copy(&graph_, &other.graph_));
  return *this;
}

inline Graph &Graph::operator=(Graph &&other) {
  graph_ = other.graph_;
#ifndef NDEBUG
  std::memset(&other.graph_, 0, sizeof(other.graph_));
#endif
  other.graph_.from.stor_begin = NULL;
  return *this;
}

inline int Graph::vcount() const noexcept { return igraph_vcount(&graph_); }

inline int Graph::ecount() const noexcept { return igraph_ecount(&graph_); }

inline bool Graph::is_directed() const noexcept {
  return igraph_is_directed(&graph_);
}

inline bool Graph::is_connected(Connectedness mode) const {
  igraph_bool_t connected;
  igraph_connectedness_t m = static_cast<igraph_connectedness_t>(mode);
  int ret = igraph_is_connected(&graph_, &connected, m);
  SafeCall(ret);
  return static_cast<bool>(connected);
}

inline int Graph::diameter(Directedness directed, bool unconnected) const {
  int length = 0;
  int ret = igraph_diameter(&graph_, &length, NULL, NULL, NULL, directed,
                            unconnected);
  SafeCall(ret);
  return length;
}

inline Vector Graph::degrees(NeighborMode mode, Loops loops) const {
  Vector v;
  bool lps = (loops == AllowLoops);
  SafeCall(igraph_degree(&graph_, v.ptr(), igraph_vss_all(),
                         static_cast<igraph_neimode_t>(mode), lps));
  return std::move(v);
}

inline Graph &Graph::AddEdge(int from, int to) {
  SafeCall(igraph_add_edge(&graph_, from, to));
  return *this;
}

inline Graph &Graph::AddEdges(const Vector &edges) {
  SafeCall(igraph_add_edges(&graph_, edges.ptr(), NULL));
  return *this;
}

inline double Graph::AveragePathLength(Directedness directed,
                                       bool unconnected) const {
  double ret;
  SafeCall(igraph_average_path_length(&graph_, &ret, directed, unconnected));
  return ret;
}

inline Graph Graph::Lattice(const Vector &dimension, int nei, Directedness dir,
                            Mutuality mutual, Periodicity periodicity) {
  igraph_t graph;
  int ret =
      igraph_lattice(&graph, dimension.ptr(), nei, dir, mutual, periodicity);
  SafeCall(ret);
  return Graph(graph);
}

inline Graph Graph::ErdosRenyiGame(int vertices, double prob, Directedness dir,
                                   Loops loops) {
  igraph_t graph;
  int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, vertices,
                                    prob, dir, loops);
  SafeCall(ret);
  return Graph(graph);
}

inline Graph Graph::ErdosRenyiGame(int vertices, int edges, Directedness dir,
                                   Loops loops) {
  igraph_t graph;
  int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, vertices,
                                    edges, dir, loops);
  SafeCall(ret);
  return Graph(graph);
}

inline Graph::Graph(const igraph_t &graph) {
  SafeCall(igraph_copy(&graph_, &graph));
}

} // namespace igraph

#endif // IGRAPHPP_GRAPH_IMPL_HPP_
