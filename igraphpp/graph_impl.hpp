#ifndef IGRAPH_GRAPH_IMPL_HPP_
#define IGRAPH_GRAPH_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

Graph::Graph(const igraph_t &graph) { SafeCall(igraph_copy(&graph_, &graph)); }

Graph::Graph(int vertices, Directedness dir) {
  SafeCall(igraph_empty(&graph_, vertices, dir));
}

Graph::Graph(const VectorView &vector, int vertices, bool directed) {
  SafeCall(igraph_create(&graph_, vector.ptr(), vertices, directed));
}

Graph::~Graph() { SafeCall(igraph_destroy(&graph_)); }

Graph::Graph(const Graph &other) {
  SafeCall(igraph_copy(&graph_, &other.graph_));
}

Graph::Graph(Graph &&other) { graph_ = std::move(other.graph_); }

Graph &Graph::operator=(const Graph &other) {
  if (this == &other)
    return *this;
  SafeCall(igraph_destroy(&graph_));
  SafeCall(igraph_copy(&graph_, &other.graph_));
  return *this;
}

Graph &Graph::operator=(const Graph &&other) {
  graph_ = std::move(other.graph_);
  return *this;
}

int Graph::vcount() const noexcept { return igraph_vcount(&graph_); }

int Graph::ecount() const noexcept { return igraph_ecount(&graph_); }

bool Graph::is_directed() const noexcept { return igraph_is_directed(&graph_); }

bool Graph::is_connected(Connectedness mode) const {
  igraph_bool_t connected;
  igraph_connectedness_t m = static_cast<igraph_connectedness_t>(mode);
  int ret = igraph_is_connected(&graph_, &connected, m);
  SafeCall(ret);
  return static_cast<bool>(connected);
}

int Graph::diameter(Directedness directed, bool unconnected) const {
  int length = 0;
  int ret = igraph_diameter(&graph_, &length, NULL, NULL, NULL, directed,
                            unconnected);
  SafeCall(ret);
  return length;
}

Vector Graph::degrees(DegreeMode mode, Loops loops) const {
  Vector v;
  bool lps = (loops == AllowLoops);
  SafeCall(igraph_degree(&graph_, v.ptr(), igraph_vss_all(), static_cast<igraph_neimode_t>(mode), lps));
  return std::move(v);
}

Graph &Graph::AddEdge(int from, int to) {
  SafeCall(igraph_add_edge(&graph_, from, to));
  return *this;
}

Graph &Graph::AddEdges(const Vector &edges) {
  SafeCall(igraph_add_edges(&graph_, edges.ptr(), NULL));
  return *this;
}

double Graph::AveragePathLength(Directedness directed, bool unconnected) const {
  double ret;
  SafeCall(igraph_average_path_length(&graph_, &ret, directed, unconnected));
  return ret;
}

Graph Graph::Lattice(const Vector &dimension, int nei, Directedness dir,
                     Mutuality mutual, Periodicity periodicity) {
  igraph_t graph;
  int ret =
      igraph_lattice(&graph, dimension.ptr(), nei, dir, mutual, periodicity);
  SafeCall(ret);
  return Graph(graph);
}

Graph Graph::ErdosRenyiGame(int vertices, double prob, Directedness dir,
                            Loops loops) {
  igraph_t graph;
  int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, vertices,
                                    prob, dir, loops);
  SafeCall(ret);
  return Graph(graph);
}

Graph Graph::ErdosRenyiGame(int vertices, int edges, Directedness dir,
                            Loops loops) {
  igraph_t graph;
  int ret = igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, vertices,
                                    edges, dir, loops);
  SafeCall(ret);
  return Graph(graph);
}

} // namespace igraph

#endif // IGRAPH_GRAPH_IMPL_HPP_
