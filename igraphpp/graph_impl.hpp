#ifndef IGRAPHPP_GRAPH_IMPL_HPP_
#define IGRAPHPP_GRAPH_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <utility>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline Graph::~Graph() {
  if (owner())
    SafeCall(igraph_destroy(ptr()));
}
inline Graph::Graph(long int vertices, Directedness dir) {
  SafeCall(igraph_empty(ptr(), vertices, dir));
}
inline Graph::Graph(const VectorView &edges, long int vertices,
                    Directedness dir) {
  SafeCall(igraph_create(ptr(), edges.ptr(), vertices, dir));
}
inline Graph::Graph(std::initializer_list<double> edges, long int vertices,
                    Directedness dir) {
  Vector vector(edges);
  SafeCall(igraph_create(ptr(), vector.ptr(), vertices, dir));
}
template <typename Iterator, typename>
Graph::Graph(Iterator edges_begin, Iterator edges_end, long int vertices,
             Directedness dir) {
  Vector vector(edges_begin, edges_end);
  SafeCall(igraph_create(ptr(), vector.ptr(), vertices, dir));
}
inline Graph::Graph(const Graph &other) {
  SafeCall(igraph_copy(&graph_, &other.graph_));
}
inline Graph::Graph(Graph &&other) {
  *ptr() = *other.ptr();
  other.disown();
}

inline Graph &Graph::operator=(const Graph &other) {
  if (this == &other)
    return *this;
  if (owner())
    SafeCall(igraph_destroy(&graph_));
  SafeCall(igraph_copy(&graph_, &other.graph_));
  return *this;
}
inline Graph &Graph::operator=(Graph &&other) {
  *ptr() = *other.ptr();
  other.disown();
  return *this;
}

inline int Graph::vcount() const noexcept { return igraph_vcount(&graph_); }
inline int Graph::ecount() const noexcept { return igraph_ecount(&graph_); }
inline Edge Graph::edge(int eid) const {
  int from, to;
  SafeCall(igraph_edge(ptr(), eid, &from, &to));
  return std::make_pair(from, to);
}
inline int Graph::eid(int from, int to, Directedness dir) const {
  int eid;
  SafeCall(igraph_get_eid(ptr(), &eid, from, to, dir, false));
  return eid;
}
// Edge eids between pairs of vertices
inline Vector Graph::pairs_eids(const VectorView &pairs, Directedness dir,
                                MultiEdges multiedges) const {
  Vector eids;
  if (multiedges == IgnoreMultiEdges) {
    SafeCall(igraph_get_eids(ptr(), eids.ptr(), pairs.ptr(), NULL, dir, false));
  } else {
    SafeCall(igraph_get_eids_multi(ptr(), eids.ptr(), pairs.ptr(), NULL, dir,
                                   false));
  }
  return eids;
}
inline Vector Graph::pairs_eids(std::initializer_list<double> pairs,
                                Directedness dir, MultiEdges multiedges) const {
  Vector vector(pairs);
  return pairs_eids(vector, dir, multiedges);
}
// Edges composing a path
inline Vector Graph::path_eids(const VectorView &path, Directedness dir,
                               MultiEdges multiedges) const {
  Vector eids;
  if (multiedges == IgnoreMultiEdges) {
    SafeCall(igraph_get_eids(ptr(), eids.ptr(), NULL, path.ptr(), dir, false));
  } else {
    SafeCall(
        igraph_get_eids_multi(ptr(), eids.ptr(), NULL, path.ptr(), dir, false));
  }
  return eids;
}
inline Vector Graph::path_eids(std::initializer_list<double> path,
                               Directedness dir, MultiEdges multiedges) const {
  Vector vector(path);
  return path_eids(vector, dir, multiedges);
}
// Concatenation of |pairs_eids()| and |path_eids()|
inline Vector Graph::eids(const VectorView &pairs, const VectorView &path,
                          Directedness dir, MultiEdges multiedges) const {
  Vector eids;
  if (multiedges == IgnoreMultiEdges) {
    SafeCall(igraph_get_eids(ptr(), eids.ptr(), pairs.ptr(), path.ptr(), dir,
                             false));
  } else {
    SafeCall(igraph_get_eids_multi(ptr(), eids.ptr(), pairs.ptr(), path.ptr(),
                                   dir, false));
  }
  return eids;
}
inline Vector Graph::eids(std::initializer_list<double> pairs,
                          std::initializer_list<double> path, Directedness dir,
                          MultiEdges multiedges) const {
  Vector vpairs(pairs);
  Vector vpath(path);
  return eids(vpairs, vpath, dir, multiedges);
}
inline Vector Graph::neighbors(int vertex, NeighborMode mode) const {
  Vector neighbors;
  SafeCall(igraph_neighbors(ptr(), neighbors.ptr(), vertex,
                            static_cast<igraph_neimode_t>(mode)));
  return neighbors;
}
inline Vector Graph::incident(int vertex, NeighborMode mode) const {
  Vector incident;
  SafeCall(igraph_incident(ptr(), incident.ptr(), vertex,
                           static_cast<igraph_neimode_t>(mode)));
  return incident;
}
inline bool Graph::is_directed() const noexcept {
  return igraph_is_directed(&graph_);
}
inline Vector Graph::degree(const VertexSelector &vids, NeighborMode mode,
                            Loops loops) const {
  Vector degrees;
  SafeCall(igraph_degree(ptr(), degrees.ptr(), *vids.ptr(),
                         static_cast<igraph_neimode_t>(mode), loops));
  return degrees;
}
inline Vector Graph::degree(std::initializer_list<double> vids,
                            NeighborMode mode, Loops loops) const {
  VertexSelector sel(vids);
  Vector degrees;
  SafeCall(igraph_degree(ptr(), degrees.ptr(), *sel.ptr(),
                         static_cast<igraph_neimode_t>(mode), loops));
  return degrees;
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

// inline Vector Graph::degrees(NeighborMode mode, Loops loops) const {
//   Vector v;
//   bool lps = (loops == AllowLoops);
//   SafeCall(igraph_degree(&graph_, v.ptr(), igraph_vss_all(),
//                          static_cast<igraph_neimode_t>(mode), lps));
//   return std::move(v);
// }

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
