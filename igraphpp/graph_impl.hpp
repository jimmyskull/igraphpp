#ifndef IGRAPHPP_GRAPH_IMPL_HPP_
#define IGRAPHPP_GRAPH_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <stdexcept>
#include <utility>

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

inline Graph::~Graph() {
  if (owner()) SafeCall(igraph_destroy(ptr()));
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
  if (this == &other) return *this;
  if (owner()) SafeCall(igraph_destroy(&graph_));
  SafeCall(igraph_copy(&graph_, &other.graph_));
  return *this;
}
inline Graph &Graph::operator=(Graph &&other) {
  if (owner()) SafeCall(igraph_destroy(&graph_));
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
inline double Graph::degree(int vertex, NeighborMode mode, Loops loops) const {
  return degree(VertexSelector::Single(vertex), mode, loops).at(0);
}

inline Graph &Graph::add_edge(int from, int to) {
  SafeCall(igraph_add_edge(&graph_, from, to));
  return *this;
}
inline Graph &Graph::add_edges(const VectorView &edges) {
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
inline Graph Graph::Lattice(const VectorView &dimension, int nei,
                            Directedness dir, Mutuality mutual,
                            Periodicity periodicity) {
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
inline Graph Graph::Barabasi(int vertices, double power,
                             const VectorView &outseq, bool outpref, double A,
                             Directedness dir, BarabasiAlgorithm algo,
                             const Graph &start_graph) {
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
inline Graph Graph::StaticFitness(int vertices, const VectorView &fitness,
                                  Loops loops, bool multiple) {
  igraph_t graph;
  SafeCall(igraph_static_fitness_game(
      &graph, vertices, const_cast<igraph_vector_t *>(fitness.ptr()), NULL,
      loops, multiple));
  return Graph(graph);
}
inline Graph Graph::StaticFitness(int vertices, const VectorView &fitness_out,
                                  const VectorView &fitness_in, Loops loops,
                                  bool multiple) {
  igraph_t graph;
  SafeCall(igraph_static_fitness_game(
      &graph, vertices, const_cast<igraph_vector_t *>(fitness_out.ptr()),
      const_cast<igraph_vector_t *>(fitness_in.ptr()), loops, multiple));
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
inline Graph Graph::CallawayTraits(int vertices, int types, int edges_per_step,
                                   const VectorView &type_dist,
                                   const Matrix &pref_matrix,
                                   Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_callaway_traits_game(
      &graph, vertices, types, edges_per_step,
      const_cast<igraph_vector_t *>(type_dist.ptr()),
      const_cast<igraph_matrix_t *>(pref_matrix.ptr()), dir));
  return Graph(graph);
}
inline Graph Graph::Establishment(int vertices, int types, int trials_per_step,
                                  const VectorView &type_dist,
                                  const Matrix &pref_matrix, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_establishment_game(
      &graph, vertices, types, trials_per_step,
      const_cast<igraph_vector_t *>(type_dist.ptr()),
      const_cast<igraph_matrix_t *>(pref_matrix.ptr()), dir));
  return Graph(graph);
}
inline Graph Graph::Preference(int vertices, int types,
                               const VectorView &type_dist, bool fixed_sizes,
                               const Matrix &pref_matrix, Directedness dir,
                               Loops loops) {
  igraph_t graph;
  SafeCall(igraph_preference_game(
      &graph, vertices, types, const_cast<igraph_vector_t *>(type_dist.ptr()),
      fixed_sizes, const_cast<igraph_matrix_t *>(pref_matrix.ptr()), NULL, dir,
      loops));
  return Graph(graph);
}
inline Graph Graph::AsymmetricPreference(int vertices, int types,
                                         const Matrix &type_dist_matrix,
                                         const Matrix &pref_matrix,
                                         Loops loops) {
  igraph_t graph;
  SafeCall(igraph_asymmetric_preference_game(
      &graph, vertices, types,
      const_cast<igraph_matrix_t *>(type_dist_matrix.ptr()),
      const_cast<igraph_matrix_t *>(pref_matrix.ptr()), NULL, NULL, loops));
  return Graph(graph);
}
inline Graph Graph::RecentDegree(int vertices, double power, int window,
                                 bool outpref, double zero_appeal,
                                 Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_recent_degree_game(&graph, vertices, power, window, 0, NULL,
                                     outpref, zero_appeal, dir));
  return Graph(graph);
}
inline Graph Graph::RecentDegree(int vertices, double power, int window, int m,
                                 const VectorView &outseq, bool outpref,
                                 double zero_appeal, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_recent_degree_game(&graph, vertices, power, window, m,
                                     outseq.ptr(), outpref, zero_appeal, dir));
  return Graph(graph);
}
inline Graph Graph::BarabasiAging(int vertices, int m, bool outpref,
                                  double pa_exp, double aging_exp,
                                  int aging_bin, double zero_deg_appeal,
                                  double zero_age_appeal, double deg_coef,
                                  double age_coef, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_barabasi_aging_game(
      &graph, vertices, m, NULL, outpref, pa_exp, aging_exp, aging_bin,
      zero_deg_appeal, zero_age_appeal, deg_coef, age_coef, dir));
  return Graph(graph);
}
inline Graph Graph::BarabasiAging(int vertices, const VectorView &outseq,
                                  bool outpref, double pa_exp, double aging_exp,
                                  int aging_bin, double zero_deg_appeal,
                                  double zero_age_appeal, double deg_coef,
                                  double age_coef, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_barabasi_aging_game(
      &graph, vertices, 0, outseq.ptr(), outpref, pa_exp, aging_exp, aging_bin,
      zero_deg_appeal, zero_age_appeal, deg_coef, age_coef, dir));
  return Graph(graph);
}
inline Graph Graph::RecentDegreeAging(int vertices, int m, bool outpref,
                                      double pa_exp, double aging_exp,
                                      int aging_bin, double time_window,
                                      double zero_appeal, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_recent_degree_aging_game(&graph, vertices, m, NULL, outpref,
                                           pa_exp, aging_exp, aging_bin,
                                           time_window, zero_appeal, dir));
  return Graph(graph);
}
inline Graph Graph::RecentDegreeAging(int vertices, const VectorView &outseq,
                                      bool outpref, double pa_exp,
                                      double aging_exp, int aging_bin,
                                      double time_window, double zero_appeal,
                                      Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_recent_degree_aging_game(
      &graph, vertices, 0, outseq.ptr(), outpref, pa_exp, aging_exp, aging_bin,
      time_window, zero_appeal, dir));
  return Graph(graph);
}
inline Graph Graph::CitedType(int vertices, const VectorView &types,
                              const VectorView &pref, int edges_per_step,
                              Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_cited_type_game(&graph, vertices, types.ptr(), pref.ptr(),
                                  edges_per_step, dir));
  return Graph(graph);
}

/* Graph, vertex, and edge attributes */

/* Structural properties of graphs */
/* Basic properties */
inline bool Graph::are_connected(int v1, int v2) const {
  int connected;
  SafeCall(igraph_are_connected(ptr(), v1, v2, &connected));
  return connected != 0;
}

/* Shortest path related functions */
inline Matrix Graph::shortest_paths(const VertexSelector &from,
                                    const VertexSelector &to,
                                    NeighborMode mode) const {
  Matrix result(from.size(*this), to.size(*this));
  SafeCall(igraph_shortest_paths(ptr(), result.ptr(), from.vs(), to.vs(),
                                 static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline double Graph::shortest_paths(int from, int to, NeighborMode mode) const {
  return shortest_paths(VertexSelector::Single(from),
                        VertexSelector::Single(to), mode)
      .at(0, 0);
}
inline Matrix Graph::shortest_paths_dijkstra(const VertexSelector &from,
                                             const VertexSelector &to,
                                             const VectorView &weights,
                                             NeighborMode mode) const {
  Matrix result(from.size(*this), to.size(*this));
  SafeCall(igraph_shortest_paths_dijkstra(ptr(), result.ptr(), from.vs(),
                                          to.vs(), weights.ptr(),
                                          static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline double Graph::shortest_paths_dijkstra(int from, int to,
                                             const VectorView &weights,
                                             NeighborMode mode) const {
  return shortest_paths_dijkstra(VertexSelector::Single(from),
                                 VertexSelector::Single(to), weights, mode)
      .at(0, 0);
}
inline Matrix Graph::shortest_paths_bellman_ford(const VertexSelector &from,
                                                 const VertexSelector &to,
                                                 const VectorView &weights,
                                                 NeighborMode mode) const {
  Matrix result(from.size(*this), to.size(*this));
  SafeCall(igraph_shortest_paths_bellman_ford(
      ptr(), result.ptr(), from.vs(), to.vs(), weights.ptr(),
      static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline double Graph::shortest_paths_bellman_ford(int from, int to,
                                                 const VectorView &weights,
                                                 NeighborMode mode) const {
  return shortest_paths_bellman_ford(VertexSelector::Single(from),
                                     VertexSelector::Single(to), weights, mode)
      .at(0, 0);
}
inline Matrix Graph::shortest_paths_johnson(const VertexSelector &from,
                                            const VertexSelector &to,
                                            const VectorView &weights) const {
  Matrix result(from.size(*this), to.size(*this));
  SafeCall(igraph_shortest_paths_johnson(ptr(), result.ptr(), from.vs(),
                                         to.vs(), weights.ptr()));
  return result;
}
inline double Graph::shortest_paths_johnson(int from, int to,
                                            const VectorView &weights) const {
  return shortest_paths_johnson(VertexSelector::Single(from),
                                VertexSelector::Single(to), weights)
      .at(0, 0);
}
// Skipped igraph_get_shortest_paths
// Skipped igraph_get_shortest_path
// Skipped igraph_get_shortest_paths_dijkstra
// Skipped igraph_get_shortest_path_dijkstra
// Skipped igraph_get_all_shortest_paths
// Skipped igraph_get_all_shortest_paths_dijkstra
inline double Graph::average_path_length(Directedness dir,
                                         bool unconnected) const {
  double ret;
  SafeCall(igraph_average_path_length(&graph_, &ret, dir, unconnected));
  return ret;
}
inline Vector Graph::path_length_hist(Directedness dir,
                                      double *unconnected) const {
  Vector result;
  SafeCall(igraph_path_length_hist(ptr(), result.ptr(), unconnected, dir));
  return result;
}
inline int Graph::diameter(int *source_vertex, int *target_vertex,
                           Directedness dir, bool unconnected) const {
  int length = 0;
  SafeCall(igraph_diameter(&graph_, &length, source_vertex, target_vertex, NULL,
                           dir, unconnected));
  return length;
}
inline int Graph::diameter(VectorView &path, Directedness dir,
                           bool unconnected) const {
  int length = 0;
  SafeCall(igraph_diameter(&graph_, &length, NULL, NULL, path.ptr(), dir,
                           unconnected));
  return length;
}
inline double Graph::diameter_dijkstra(const VectorView &weights,
                                       int *source_vertex, int *target_vertex,
                                       Directedness dir,
                                       bool unconnected) const {
  double result = 0;
  SafeCall(igraph_diameter_dijkstra(&graph_, weights.ptr(), &result,
                                    source_vertex, target_vertex, NULL, dir,
                                    unconnected));
  return result;
}
inline double Graph::diameter_dijkstra(const VectorView &weights,
                                       VectorView &path, Directedness dir,
                                       bool unconnected) const {
  double result = 0;
  SafeCall(igraph_diameter_dijkstra(&graph_, weights.ptr(), &result, NULL, NULL,
                                    path.ptr(), dir, unconnected));
  return result;
}
inline int Graph::girth() const {
  int girth;
  SafeCall(igraph_girth(ptr(), &girth, NULL));
  return girth;
}
inline int Graph::girth(VectorView &circle) const {
  int girth;
  SafeCall(igraph_girth(ptr(), &girth, circle.ptr()));
  return girth;
}
inline Vector Graph::eccentricity(const VertexSelector &vids,
                                  NeighborMode mode) const {
  Vector result;
  SafeCall(igraph_eccentricity(ptr(), result.ptr(), vids.vs(),
                               static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline double Graph::eccentricity(int vertex, NeighborMode mode) const {
  return eccentricity(VertexSelector::Single(vertex), mode).at(0);
}
inline double Graph::radius(NeighborMode mode) const {
  double radius;
  SafeCall(igraph_radius(ptr(), &radius, static_cast<igraph_neimode_t>(mode)));
  return radius;
}

/* Neighborhood of a vertex */
inline Vector Graph::neighborhood_size(const VertexSelector &vids, int order,
                                       NeighborMode mode) const {
  Vector result;
  SafeCall(igraph_neighborhood_size(ptr(), result.ptr(), vids.vs(), order,
                                    static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline double Graph::neighborhood_size(int vertex, int order,
                                       NeighborMode mode) const {
  return neighborhood_size(VertexSelector::Single(vertex), order, mode).at(0);
}

/* Graph components */
inline Vector Graph::subcomponent(int vertex, NeighborMode mode) const {
  Vector result;
  SafeCall(igraph_subcomponent(ptr(), result.ptr(), static_cast<double>(vertex),
                               static_cast<igraph_neimode_t>(mode)));
  return result;
}
inline Graph Graph::induced_subgraph(
    const VertexSelector &vids, SubgraphImplementation implementation) const {
  igraph_t graph;
  SafeCall(igraph_induced_subgraph(
      ptr(), &graph, vids.vs(),
      static_cast<igraph_subgraph_implementation_t>(implementation)));
  return Graph(graph);
}
inline Graph Graph::subgraph_edges(const EdgeSelector &eids,
                                   bool delete_vertices) const {
  igraph_t graph;
  SafeCall(igraph_subgraph_edges(ptr(), &graph, eids.es(), delete_vertices));
  return Graph(graph);
}
inline Clusters Graph::clusters(Connectedness mode) const {
  Clusters clusters;
  SafeCall(igraph_clusters(ptr(), clusters.membership.ptr(),
                           clusters.csize.ptr(), &clusters.number_of_clusters,
                           static_cast<igraph_connectedness_t>(mode)));
  return clusters;
}
inline bool Graph::is_connected(Connectedness mode) const {
  int connected;
  SafeCall(igraph_is_connected(&graph_, &connected,
                               static_cast<igraph_connectedness_t>(mode)));
  return static_cast<bool>(connected);
}
inline Vector Graph::articulation_points() const {
  Vector vector;
  SafeCall(igraph_articulation_points(ptr(), vector.ptr()));
  return vector;
}

/* Degree sequences */
inline bool Graph::is_degree_sequence(const Vector &out_degrees,
                                      const Vector &in_degrees) {
  int result;
  SafeCall(
      igraph_is_degree_sequence(out_degrees.ptr(), in_degrees.ptr(), &result));
  return result != 0;
}
inline bool Graph::is_graphical_degree_sequence(const Vector &out_degrees,
                                                const Vector &in_degrees) {
  int result;
  SafeCall(igraph_is_graphical_degree_sequence(out_degrees.ptr(),
                                               in_degrees.ptr(), &result));
  return result != 0;
}

/* Centrality measures */
inline Vector Graph::closeness(const VertexSelector &vids, NeighborMode mode,
                               const VectorView &weights,
                               bool normalized) const {
  Vector result;
  SafeCall(igraph_closeness(
      ptr(), result.ptr(), vids.vs(), static_cast<igraph_neimode_t>(mode),
      weights.is_none() ? NULL : weights.ptr(), normalized));
  return result;
}
inline double Graph::closeness(int vertex, NeighborMode mode,
                               const VectorView &weights,
                               bool normalized) const {
  return closeness(VertexSelector::Single(vertex), mode, weights, normalized)
      .at(0);
}
inline Vector Graph::betweenness(const VertexSelector &vids, Directedness dir,
                                 const VectorView &weights,
                                 bool nobigint) const {
  Vector result;
  SafeCall(igraph_betweenness(ptr(), result.ptr(), vids.vs(), dir,
                              weights.ptr(), nobigint));
  return result;
}
inline double Graph::betweenness(int vertex, Directedness dir,
                                 const VectorView &weights,
                                 bool nobigint) const {
  return betweenness(VertexSelector::Single(vertex), dir, weights, nobigint)
      .at(0);
}
inline Vector Graph::edge_betweenness(Directedness dir,
                                      const VectorView &weights) const {
  Vector result;
  SafeCall(igraph_edge_betweenness(ptr(), result.ptr(), dir, weights.ptr()));
  return result;
}
inline PageRank Graph::pagerank(const VertexSelector &vids, Directedness dir,
                                double damping,
                                const VectorView &weights) const {
  PageRank result;
  SafeCall(igraph_pagerank(ptr(), IGRAPH_PAGERANK_ALGO_PRPACK,
                           result.scores.ptr(), &result.eigenvalue, vids.vs(),
                           dir, damping, weights.ptr(), NULL));
  return result;
}
inline Vector Graph::constraint(const VertexSelector &vids,
                                const VectorView &weights) const {
  Vector result;
  SafeCall(igraph_constraint(ptr(), result.ptr(), vids.vs(), weights.ptr()));
  return result;
}
inline double Graph::constraint(int vertex, const VectorView &weights) const {
  return constraint(VertexSelector::Single(vertex), weights).at(0);
}
inline int Graph::maxdegree(const VertexSelector &vids, NeighborMode mode,
                            Loops loops) const {
  int result;
  SafeCall(igraph_maxdegree(ptr(), &result, vids.vs(),
                            static_cast<igraph_neimode_t>(mode), loops));
  return result;
}
inline Vector Graph::strength(const VertexSelector &vids, NeighborMode mode,
                              Loops loops, const VectorView &weights) const {
  Vector result;
  SafeCall(igraph_strength(ptr(), result.ptr(), vids.vs(),
                           static_cast<igraph_neimode_t>(mode), loops,
                           weights.ptr()));
  return result;
}
inline double Graph::strength(int vertex, NeighborMode mode, Loops loops,
                              const VectorView &weights) const {
  return strength(VertexSelector::Single(vertex), mode, loops, weights).at(0);
}

inline Graph Graph::ReadEdgelist(FILE *instream, int n, Directedness dir) {
  igraph_t graph;
  SafeCall(igraph_read_graph_edgelist(&graph, instream, n, dir));
  return Graph(graph);
}
inline Graph Graph::ReadEdgelist(std::string filename, int n, Directedness dir) {
  FILE *fp = fopen(filename.c_str(), "r");
  if (fp == NULL)
    throw std::runtime_error("File not found");
  Graph graph = ReadEdgelist(fp, n, dir);
  fclose(fp);
  return graph;
}


inline Graph::Graph(const igraph_t &graph) : graph_(graph) {}

}  // namespace igraph

#endif  // IGRAPHPP_GRAPH_IMPL_HPP_
