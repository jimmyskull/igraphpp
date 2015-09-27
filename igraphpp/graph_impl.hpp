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
  if (owner())
    SafeCall(igraph_destroy(&graph_));
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
  SafeCall(igraph_degree(ptr(), degrees.ptr(), vids.vs(),
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

inline Graph &Graph::add_edge(int from, int to) {
  SafeCall(igraph_add_edge(&graph_, from, to));
  return *this;
}
inline Graph &Graph::add_edges(const Vector &edges) {
  SafeCall(igraph_add_edges(&graph_, edges.ptr(), NULL));
  return *this;
}
inline Graph &Graph::add_edges(std::initializer_list<double> edges) {
  Vector vector(edges);
  return add_edges(vector);
}
inline Graph &Graph::add_vertices(int number_of_vertices) {
  SafeCall(igraph_add_vertices(ptr(), number_of_vertices, NULL));
  return *this;
}
inline void Graph::delete_edges(const EdgeSelector &edges) {
  SafeCall(igraph_delete_edges(ptr(), edges.es()));
}
inline void Graph::delete_edges(std::initializer_list<double> edges) {
  EdgeSelector sel(edges);
  SafeCall(igraph_delete_edges(ptr(), sel.es()));
}
inline void Graph::delete_vertices(const VertexSelector &vertices) {
  SafeCall(igraph_delete_vertices(ptr(), vertices.vs()));
}
inline void Graph::delete_vertices(std::initializer_list<double> vertices) {
  VertexSelector sel(vertices);
  SafeCall(igraph_delete_vertices(ptr(), sel.vs()));
}

inline int Graph::diameter(Directedness directed, bool unconnected) const {
  int length = 0;
  int ret = igraph_diameter(&graph_, &length, NULL, NULL, NULL, directed,
                            unconnected);
  SafeCall(ret);
  return length;
}

inline bool Graph::is_connected(Connectedness mode) const {
  igraph_bool_t connected;
  igraph_connectedness_t m = static_cast<igraph_connectedness_t>(mode);
  int ret = igraph_is_connected(&graph_, &connected, m);
  SafeCall(ret);
  return static_cast<bool>(connected);
}

inline double Graph::AveragePathLength(Directedness directed,
                                       bool unconnected) const {
  double ret;
  SafeCall(igraph_average_path_length(&graph_, &ret, directed, unconnected));
  return ret;
}

inline Graph Graph::AdjacencyMatrix(Matrix &adjmatrix,
                                    AdjacencyMatrixMode mode) {
  igraph_t graph;
  SafeCall(igraph_adjacency(&graph, adjmatrix.ptr(),
                            static_cast<igraph_adjacency_t>(mode)));
  return Graph(graph);
}
inline Graph Graph::Star(int vertices, StarMode mode, int center_vertex) {
  igraph_t graph;
  SafeCall(igraph_star(&graph, vertices, static_cast<igraph_star_mode_t>(mode),
                       center_vertex));
  return Graph(graph);
}
inline Graph Graph::Lattice(const Vector &dimension, int nei, Directedness dir,
                            Mutuality mutual, Periodicity periodicity) {
  igraph_t graph;
  SafeCall(
      igraph_lattice(&graph, dimension.ptr(), nei, dir, mutual, periodicity));
  return Graph(graph);
}
inline Graph Graph::Ring(int vertices, Directedness dir, Mutuality mutual,
                         bool circular) {
  igraph_t graph;
  SafeCall(igraph_ring(&graph, vertices, dir, mutual, circular));
  return Graph(graph);
}
inline Graph Graph::Tree(int vertices, int children, TreeMode mode) {
  igraph_t graph;
  SafeCall(igraph_tree(&graph, vertices, children,
                       static_cast<igraph_tree_mode_t>(mode)));
  return Graph(graph);
}
inline Graph Graph::Full(int vertices, Directedness dir, Loops loops) {
  igraph_t graph;
  SafeCall(igraph_full(&graph, vertices, dir, loops));
  return Graph(graph);
}
inline Graph Graph::FullCitation(int vertices, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_full_citation(&graph, vertices, dir));
  return Graph(graph);
}
inline Graph Graph::Famous(const char *name) {
  igraph_t graph;
  SafeCall(igraph_famous(&graph, name));
  return Graph(graph);
}
template <typename... Args, typename>
Graph Graph::LCF(int vertices, Args... args) {
  igraph_t graph;
  SafeCall(igraph_lcf(&graph, vertices, args..., 0));
  return Graph(graph);
}
inline Graph Graph::LCF(int vertices, const VectorView &shifts, int repeats) {
  igraph_t graph;
  SafeCall(igraph_lcf_vector(&graph, vertices, shifts.ptr(), repeats));
  return Graph(graph);
}
inline Graph Graph::LCF(int vertices, std::initializer_list<double> shifts,
                        int repeats) {
  Vector vshifts(shifts);
  return LCF(vertices, vshifts, repeats);
}
inline Graph Graph::Atlas(int number) {
  igraph_t graph;
  SafeCall(igraph_atlas(&graph, number));
  return Graph(graph);
}
inline Graph Graph::deBruijn(int m, int n) {
  igraph_t graph;
  SafeCall(igraph_de_bruijn(&graph, m, n));
  return Graph(graph);
}
inline Graph Graph::Kautz(int m, int n) {
  igraph_t graph;
  SafeCall(igraph_kautz(&graph, m, n));
  return Graph(graph);
}
inline Graph Graph::ExtendedChordalRing(int vertices, const Matrix &W) {
  igraph_t graph;
  SafeCall(igraph_extended_chordal_ring(&graph, vertices, W.ptr()));
  return Graph(graph);
}
inline void Graph::connect_neighborhood(int order, NeighborMode mode) {
  SafeCall(igraph_connect_neighborhood(ptr(), order,
                                       static_cast<igraph_neimode_t>(mode)));
}

/* Randomized graph generators */
inline Graph Graph::GRG(int vertices, double radius, bool torus) {
  igraph_t graph;
  SafeCall(igraph_grg_game(&graph, vertices, radius, torus, NULL, NULL));
  return Graph(graph);
}
inline Graph Graph::Barabasi(int vertices, double power, int m, bool outpref,
                             double A, Directedness dir, BarabasiAlgorithm algo,
                             const Graph &start_graph) {
  igraph_t graph;
  const igraph_t *start = start_graph.vcount() ? start_graph.ptr() : NULL;
  SafeCall(igraph_barabasi_game(
      &graph, vertices, power, m, NULL, outpref, A, dir,
      static_cast<igraph_barabasi_algorithm_t>(algo), start));
  return Graph(graph);
}
inline Graph Graph::Barabasi(int vertices, double power, const Vector &outseq,
                             bool outpref, double A, Directedness dir,
                             BarabasiAlgorithm algo, const Graph &start_graph) {
  igraph_t graph;
  const igraph_t *start = start_graph.vcount() ? start_graph.ptr() : NULL;
  SafeCall(igraph_barabasi_game(
      &graph, vertices, power, 0, outseq.ptr(), outpref, A, dir,
      static_cast<igraph_barabasi_algorithm_t>(algo), start));
  return Graph(graph);
}
inline Graph Graph::ErdosRenyi(int vertices, double prob, Directedness dir,
                               Loops loops) {
  igraph_t graph;
  SafeCall(igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, vertices,
                                   prob, dir, loops));
  return Graph(graph);
}
inline Graph Graph::ErdosRenyi(int vertices, int edges, Directedness dir,
                               Loops loops) {
  igraph_t graph;
  SafeCall(igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, vertices,
                                   edges, dir, loops));
  return Graph(graph);
}
inline Graph Graph::WattsStrogatz(int dim, int size, int nei, double p,
                                  Loops loops, bool multiple) {
  igraph_t graph;
  SafeCall(
      igraph_watts_strogatz_game(&graph, dim, size, nei, p, loops, multiple));
  return Graph(graph);
}
inline void Graph::rewire_edges(double prob, Loops loops, bool multiple) {
  SafeCall(igraph_rewire_edges(ptr(), prob, loops, multiple));
}
inline Graph Graph::DegreeSequence(const VectorView &degrees,
                                   DegreeSequenceMethod method) {
  igraph_t graph;
  SafeCall(igraph_degree_sequence_game(&graph, degrees.ptr(), NULL,
                                       static_cast<igraph_degseq_t>(method)));
  return Graph(graph);
}
inline Graph Graph::DegreeSequence(const VectorView &out_deq,
                                   const VectorView &in_deq,
                                   DegreeSequenceMethod method) {
  igraph_t graph;
  SafeCall(igraph_degree_sequence_game(&graph, out_deq.ptr(), in_deq.ptr(),
                                       static_cast<igraph_degseq_t>(method)));
  return Graph(graph);
}
inline Graph Graph::kRegular(int vertices, int k, Directedness dir,
                             bool multiple) {
  igraph_t graph;
  SafeCall(igraph_k_regular_game(&graph, vertices, k, dir, multiple));
  return Graph(graph);
}
inline Graph Graph::StaticFitness(int vertices, VectorView &fitness,
                                  Loops loops, bool multiple) {
  igraph_t graph;
  SafeCall(igraph_static_fitness_game(&graph, vertices, fitness.ptr(), NULL,
                                      loops, multiple));
  return Graph(graph);
}
inline Graph Graph::StaticFitness(int vertices, VectorView &fitness_out,
                                  VectorView &fitness_in, Loops loops,
                                  bool multiple) {
  igraph_t graph;
  SafeCall(igraph_static_fitness_game(&graph, vertices, fitness_out.ptr(),
                                      fitness_in.ptr(), loops, multiple));
  return Graph(graph);
}
inline Graph Graph::StaticPowerLaw(int vertices, int edges, double exponent_out,
                                   double exponent_in, Loops loops,
                                   bool multiple,
                                   bool finite_size_correlation) {
  igraph_t graph;
  SafeCall(igraph_static_power_law_game(&graph, vertices, edges, exponent_out,
                                        exponent_in, loops, multiple,
                                        finite_size_correlation));
  return Graph(graph);
}
inline Graph Graph::ForestFire(int vertices, double fw_prob, double bw_factor,
                               int pambs, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_forest_fire_game(&graph, vertices, fw_prob, bw_factor, pambs,
                                   dir));
  return Graph(graph);
}
inline void Graph::rewire(int trials, RewiringMode mode) {
  SafeCall(igraph_rewire(ptr(), trials, static_cast<igraph_rewiring_t>(mode)));
}
inline Graph Graph::GrowingRandom(int vertices, int m, Directedness dir,
                           bool citation) {
  igraph_t graph;
  SafeCall(igraph_growing_random_game(&graph, vertices, m, dir, citation));
  return Graph(graph);
}

inline Graph::Graph(const igraph_t &graph) : graph_(graph) {}

} // namespace igraph

#endif // IGRAPHPP_GRAPH_IMPL_HPP_
