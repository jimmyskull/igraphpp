#ifndef IGRAPHPP_GRAPH_HPP_
#define IGRAPHPP_GRAPH_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#ifndef NDEBUG
#include <cstring>
#endif

#include <initializer_list>
#include <functional>
#include <string>

#include <igraph/igraph.h>

#include "./exception.hpp"
#include "./mapper.hpp"
#include "./vector.hpp"
#include "./util.hpp"

namespace igraph {

class Vector;
class VertexSelector;

struct Clusters {
  Vector membership;
  Vector csize;
  int number_of_clusters;
};

struct PageRank {
  Vector scores;
  double eigenvalue;
};

auto bfs_default_callback = [](int /*vid*/, int /*pred*/, int /*succ*/,
                               int /*rank*/,
                               int /*dist*/) -> bool { return false; };

class Graph {
  template <typename>
  friend struct TypeMapper;

 public:
  /* Constructors and Destructors */
  ~Graph();
  Graph(long int vertices = 0, bool directed = Undirected);
  // TODO: Empty graph with attributes
  Graph(const VectorView &edges, long int vertices = 0,
        bool directed = Undirected);
  Graph(std::initializer_list<double> edges, long int vertices = 0,
        bool directed = Undirected);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  Graph(Iterator edges_begin, Iterator edges_end, long int vertices = 0,
        bool directed = Undirected);
  Graph(const Graph &other);
  Graph(Graph &&other);

  Graph &operator=(const Graph &other)&;
  Graph &operator=(Graph &&other)&;

  /* Basic query operations */
  int vcount() const noexcept;
  int ecount() const noexcept;
  Edge edge(int eid) const;
  int eid(int from, int to, bool directed = Undirected) const;
  // Edge eids between pairs of vertices
  Vector pairs_eids(const VectorView &pairs, bool directed = Directed,
                    MultiEdges multiedges = MultiEdges::Ignore) const;
  Vector pairs_eids(std::initializer_list<double> pairs,
                    bool directed = Directed,
                    MultiEdges multiedges = MultiEdges::Ignore) const;
  // Edges composing a path
  Vector path_eids(const VectorView &path, bool directed = Directed,
                   MultiEdges multiedges = MultiEdges::Ignore) const;
  Vector path_eids(std::initializer_list<double> path, bool directed = Directed,
                   MultiEdges multiedges = MultiEdges::Ignore) const;
  // Concatenation of |pairs_eids()| and |path_eids()|
  Vector eids(const VectorView &pairs, const VectorView &path,
              bool directed = Directed,
              MultiEdges multiedges = MultiEdges::Ignore) const;
  Vector eids(std::initializer_list<double> pairs,
              std::initializer_list<double> path, bool directed = Directed,
              MultiEdges multiedges = MultiEdges::Ignore) const;
  Vector neighbors(int vertex, Mode mode = Mode::Out) const;
  Vector incident(int vertex, Mode mode = Mode::Out) const;
  bool is_directed() const noexcept;
  Vector degree(const VertexSelector &vids = VertexSelector::All(),
                Mode mode = Mode::Out, Loops loops = Loops::None) const;
  void degree(Vector *degrees,
              const VertexSelector &vids = VertexSelector::All(),
              Mode mode = Mode::Out, Loops loops = Loops::None) const;
  double degree(int vertex, Mode mode = Mode::Out,
                Loops loops = Loops::None) const;

  /* Adding and deleting vertices and edges */
  Graph &add_edge(int from, int to);
  Graph &add_edges(const VectorView &edges);
  Graph &add_edges(std::initializer_list<double> edges);
  Graph &add_vertices(int number_of_vertices);
  void delete_edges(const EdgeSelector &edges);
  void delete_edges(std::initializer_list<double> edges);
  void delete_vertices(const VertexSelector &vertices);
  void delete_vertices(std::initializer_list<double> vertices);

  /* Deterministic graph generators */
  static Graph AdjacencyMatrix(
      Matrix &adjmatrix,
      AdjacencyMatrixMode mode = AdjacencyMatrixMode::Directed);
  // Skipped igraph_weighted_adjacency (attributes has not been implemented)
  // Skipped  igraph_adjlist (igraph_adjlist_t has not been implemented)
  static Graph Star(int vertices, StarMode mode = StarMode::Out,
                    int center_vertex = 0);
  static Graph Lattice(const VectorView &dimension, int nei = 1,
                       bool directed = Undirected,
                       Mutuality mutual = Mutuality::None,
                       Periodicity periodicity = Periodicity::None);
  static Graph Ring(int vertices, bool directed = Undirected,
                    Mutuality mutual = Mutuality::None, bool circular = true);
  static Graph Tree(int vertices, int children = 2,
                    TreeMode mode = TreeMode::Out);
  static Graph Full(int vertices, bool directed = Undirected,
                    Loops loops = Loops::None);
  static Graph FullCitation(int vertices, bool directed = Directed);
  static Graph Famous(const char *name);
  template <typename... Args, typename = typename std::enable_if<util::all_args(
                                  std::is_same<Args, int>::value...)>::type>
  static Graph LCF(int vertices, Args... args);
  static Graph LCF(int vertices, const VectorView &shifts, int repeats);
  static Graph LCF(int vertices, std::initializer_list<double> shifts,
                   int repeats);
  static Graph Atlas(int number);
  static Graph deBruijn(int m, int n);
  static Graph Kautz(int m, int n);
  static Graph ExtendedChordalRing(int vertices, const Matrix &W);
  void connect_neighborhood(int order, Mode mode = Mode::Out);

  /* Randomized graph generators */
  static Graph GRG(int vertices, double radius, bool torus);
  static Graph Barabasi(int vertices, double power = 1, int m = 2,
                        bool outpref = false, double A = 1,
                        bool directed = Directed,
                        BarabasiAlgorithm algo = BarabasiAlgorithm::PSumTree,
                        const Graph &start_graph = Graph());
  static Graph Barabasi(int vertices, double power, const VectorView &outseq,
                        bool outpref = false, double A = 1,
                        bool directed = Directed,
                        BarabasiAlgorithm algo = BarabasiAlgorithm::PSumTree,
                        const Graph &start_graph = Graph());
  static Graph ErdosRenyi(int vertices, double prob, bool directed = Undirected,
                          Loops loops = Loops::None);
  static Graph ErdosRenyi(int vertices, int edges, bool directed = Undirected,
                          Loops loops = Loops::None);
  static Graph WattsStrogatz(int dim, int size, int nei, double p,
                             Loops loops = Loops::None, bool multiple = false);
  // Rewire uniformly
  void rewire_edges(double prob, Loops loops = Loops::None,
                    bool multiple = false);
  static Graph DegreeSequence(
      const VectorView &degrees,
      DegreeSequenceMethod method = DegreeSequenceMethod::Simple);
  static Graph DegreeSequence(
      const VectorView &out_deq, const VectorView &in_deq,
      DegreeSequenceMethod method = DegreeSequenceMethod::Simple);
  static Graph kRegular(int vertices, int k, bool directed = Directed,
                        bool multiple = false);
  static Graph StaticFitness(int vertices, const VectorView &fitness,
                             Loops loops = Loops::None, bool multiple = false);
  static Graph StaticFitness(int vertices, const VectorView &fitness_out,
                             const VectorView &fitness_in,
                             Loops loops = Loops::None, bool multiple = false);
  static Graph StaticPowerLaw(int vertices, int edges, double exponent_out,
                              double exponent_in = -1.0,
                              Loops loops = Loops::None, bool multiple = false,
                              bool finite_size_correlation = true);
  static Graph ForestFire(int vertices, double fw_prob, double bw_factor = 1.0,
                          int pambs = 1, bool directed = Undirected);
  // Rewire keeping the degree distribution
  void rewire(int trials, RewiringMode mode = RewiringMode::Simple);
  static Graph GrowingRandom(int vertices, int m = 1, bool directed = Directed,
                             bool citation = false);
  static Graph CallawayTraits(int vertices, int types, int edges_per_step,
                              const VectorView &type_dist,
                              const Matrix &pref_matrix,
                              bool directed = Undirected);
  static Graph Establishment(int vertices, int types, int trials_per_step,
                             const VectorView &type_dist,
                             const Matrix &pref_matrix,
                             bool directed = Undirected);
  static Graph Preference(int vertices, int types, const VectorView &type_dist,
                          bool fixed_sizes, const Matrix &pref_matrix,
                          bool directed = Undirected,
                          Loops loops = Loops::None);
  static Graph AsymmetricPreference(int vertices, int types,
                                    const Matrix &type_dist_matrix,
                                    const Matrix &pref_matrix,
                                    Loops loops = Loops::None);
  static Graph RecentDegree(int vertices, double power, int window,
                            bool outpref = false, double zero_appeal = 1.0,
                            bool directed = Undirected);
  static Graph RecentDegree(int vertices, double power, int window, int m,
                            const VectorView &outseq, bool outpref = false,
                            double zero_appeal = 1.0,
                            bool directed = Undirected);
  static Graph BarabasiAging(int vertices, int m, bool outpref, double pa_exp,
                             double aging_exp, int aging_bin,
                             double zero_deg_appeal, double zero_age_appeal,
                             double deg_coef, double age_coef,
                             bool directed = Undirected);
  static Graph BarabasiAging(int vertices, const VectorView &outseq,
                             bool outpref, double pa_exp, double aging_exp,
                             int aging_bin, double zero_deg_appeal,
                             double zero_age_appeal, double deg_coef,
                             double age_coef, bool directed = Undirected);
  static Graph RecentDegreeAging(int vertices, int m, bool outpref,
                                 double pa_exp, double aging_exp, int aging_bin,
                                 double time_window, double zero_appeal,
                                 bool directed = Undirected);
  static Graph RecentDegreeAging(int vertices, const VectorView &outseq,
                                 bool outpref, double pa_exp, double aging_exp,
                                 int aging_bin, double time_window,
                                 double zero_appeal,
                                 bool directed = Undirected);
  static Graph CitedType(int vertices, const VectorView &types,
                         const VectorView &pref, int edges_per_step,
                         bool directed = Undirected);
  // Skipped igraph_sbm_game (igraph_vector_int_t has not been implemented)

  /* Graph, vertex, and edge attributes */

  /* Structural properties of graphs */
  /* Basic properties */
  bool are_connected(int v1, int v2) const;

  /* Shortest path related functions */
  Matrix shortest_paths(const VertexSelector &from, const VertexSelector &to,
                        Mode mode = Mode::Out) const;
  double shortest_paths(int from, int to, Mode mode = Mode::Out) const;
  Matrix shortest_paths_dijkstra(const VertexSelector &from,
                                 const VertexSelector &to,
                                 const VectorView &weights,
                                 Mode mode = Mode::Out) const;
  double shortest_paths_dijkstra(int from, int to, const VectorView &weights,
                                 Mode mode = Mode::Out) const;
  Matrix shortest_paths_bellman_ford(const VertexSelector &from,
                                     const VertexSelector &to,
                                     const VectorView &weights,
                                     Mode mode = Mode::Out) const;
  double shortest_paths_bellman_ford(int from, int to,
                                     const VectorView &weights,
                                     Mode mode = Mode::Out) const;
  Matrix shortest_paths_johnson(const VertexSelector &from,
                                const VertexSelector &to,
                                const VectorView &weights) const;
  double shortest_paths_johnson(int from, int to,
                                const VectorView &weights) const;
  // Skipped igraph_get_shortest_paths
  // Skipped igraph_get_shortest_path
  // Skipped igraph_get_shortest_paths_dijkstra
  // Skipped igraph_get_shortest_path_dijkstra
  // Skipped igraph_get_all_shortest_paths
  // Skipped igraph_get_all_shortest_paths_dijkstra
  double average_path_length(bool directed = Directed,
                             bool unconnected = true) const;
  Vector path_length_hist(bool directed = Directed,
                          double *unconnected = NULL) const;
  int diameter(int *source_vertex = NULL, int *target_vertex = NULL,
               bool directed = Directed, bool unconnected = true) const;
  int diameter(VectorView &path, bool directed = Directed,
               bool unconnected = true) const;
  double diameter_dijkstra(const VectorView &weights, int *source_vertex = NULL,
                           int *target_vertex = NULL, bool directed = Directed,
                           bool unconnected = true) const;
  double diameter_dijkstra(const VectorView &weights, VectorView &path,
                           bool directed = Directed,
                           bool unconnected = true) const;
  int girth() const;
  int girth(VectorView &circle) const;
  Vector eccentricity(const VertexSelector &vids = VertexSelector::All(),
                      Mode mode = Mode::Out) const;
  double eccentricity(int vertex, Mode mode = Mode::Out) const;
  double radius(Mode mode = Mode::Out) const;

  /* Neighborhood of a vertex */
  Vector neighborhood_size(const VertexSelector &vids = VertexSelector::All(),
                           int order = 1, Mode mode = Mode::Out) const;
  double neighborhood_size(int vertex, int order = 1,
                           Mode mode = Mode::Out) const;
  VectorPtr<VectorView> neighborhood(
      const VertexSelector &vids = VertexSelector::All(), int order = 1,
      Mode mode = Mode::Out);
  VectorPtr<Graph> neighborhood_graphs(
      const VertexSelector &vids = VertexSelector::All(), int order = 1,
      Mode mode = Mode::Out);

  /* Graph components */
  Vector subcomponent(int vertex, Mode mode = Mode::Out) const;
  Graph induced_subgraph(const VertexSelector &vids = VertexSelector::All(),
                         SubgraphImplementation implementation =
                             SubgraphImplementation::CreateFromScratch) const;
  Graph subgraph_edges(const EdgeSelector &eids,
                       bool delete_vertices = true) const;
  Clusters clusters(Connectedness mode = Connectedness::Weak) const;
  bool is_connected(Connectedness mode = Connectedness::Weak) const;
  // Skipped igraph_decompose
  // Skipped igraph_decompose_destroy
  // Skipped igraph_biconnected_components
  Vector articulation_points() const;

  /* Degree sequences */
  static bool is_degree_sequence(const Vector &out_degrees,
                                 const Vector &in_degrees);
  static bool is_graphical_degree_sequence(const Vector &out_degrees,
                                           const Vector &in_degrees);

  /* Centrality measures */
  Vector closeness(const VertexSelector &vids = VertexSelector::All(),
                   Mode mode = Mode::Out,
                   const VectorView &weights = VectorView::None(),
                   bool normalized = false) const;
  double closeness(int vertex, Mode mode = Mode::Out,
                   const VectorView &weights = VectorView::None(),
                   bool normalized = false) const;
  Vector betweenness(const VertexSelector &vids = VertexSelector::All(),
                     bool directed = Directed,
                     const VectorView &weights = VectorView::None(),
                     bool nobigint = true) const;
  double betweenness(int vertex, bool directed = Directed,
                     const VectorView &weights = VectorView::None(),
                     bool nobigint = true) const;
  Vector edge_betweenness(bool directed = Directed,
                          const VectorView &weights = VectorView::None()) const;
  PageRank pagerank(const VertexSelector &vids = VertexSelector::All(),
                    bool directed = Directed, double damping = 0.85,
                    const VectorView &weights = VectorView::None()) const;
  // Skipped igraph_personalized_pagerank
  // Skipped igraph_personalized_pagerank_vs
  Vector constraint(const VertexSelector &vids = VertexSelector::All(),
                    const VectorView &weights = VectorView::None()) const;
  double constraint(int vertex,
                    const VectorView &weights = VectorView::None()) const;
  int maxdegree(const VertexSelector &vids = VertexSelector::All(),
                Mode mode = Mode::Out, Loops loops = Loops::None) const;
  Vector strength(const VertexSelector &vids = VertexSelector::All(),
                  Mode mode = Mode::Out, Loops loops = Loops::None,
                  const VectorView &weights = VectorView::None()) const;
  double strength(int vertex, Mode mode = Mode::Out, Loops loops = Loops::None,
                  const VectorView &weights = VectorView::None()) const;
  // igraph_eigenvector_centrality
  // igraph_hub_score
  // igraph_authority_score

  /* Estimating Centrality Measures */
  // igraph_closeness_estimate
  // igraph_betweenness_estimate
  // igraph_edge_betweenness_estimate

  /* Centralization */
  // igraph_centralization
  // igraph_centralization_degree
  // igraph_centralization_betweenness
  // igraph_centralization_closeness
  // igraph_centralization_eigenvector_centrality
  // igraph_centralization_degree_tmax
  // igraph_centralization_betweenness_tmax
  // igraph_centralization_closeness_tmax
  // igraph_centralization_eigenvector_centrality_tmax

  /* Similarity Measures */
  // igraph_bibcoupling
  // igraph_cocitation
  // igraph_similarity_jaccard
  // igraph_similarity_jaccard_pairs
  // igraph_similarity_jaccard_es
  // igraph_similarity_dice
  // igraph_similarity_dice_pairs
  // igraph_similarity_dice_es
  // igraph_similarity_inverse_log_weighted

  /* Spanning Trees */
  // igraph_minimum_spanning_tree
  // igraph_minimum_spanning_tree_unweighted
  // igraph_minimum_spanning_tree_prim

  /* Transitivity or Clustering Coefficient */
  // igraph_transitivity_undirected
  // igraph_transitivity_local_undirected
  // igraph_transitivity_avglocal_undirected
  // igraph_transitivity_barrat

  /* Directedness conversion */
  Graph &to_directed(DirectedMode mode = DirectedMode::Mutual);
  Graph &to_undirected(UndirectedMode mode = UndirectedMode::Collapse);

  /* Spectral properties */
  // igraph_laplacian

  /* Non-simple graphs: multiple and loop edges */
  // igraph_is_simple
  // Skipped igraph_is_loop
  // igraph_is_multiple
  // igraph_has_multiple
  // igraph_count_multiple
  // igraph_simplify

  /* Mixing patterns */
  // igraph_assortativity_nominal
  // igraph_assortativity
  // igraph_assortativity_degree

  /* K-Cores */
  // igraph_coreness

  /* Topological sorting, directed acyclic graphs */
  // igraph_is_dag
  // igraph_topological_sorting
  // igraph_feedback_arc_set

  /* Maximum cardinality search, graph decomposition, chordal graphs */
  // igraph_maximum_cardinality_search
  // igraph_is_chordal

  /* Matchings */
  // Skipped igraph_is_matching
  // Skipped igraph_is_maximal_matching
  // Skipped igraph_maximum_bipartite_matching

  /* Line graphs */
  // igraph_linegraph

  /* Unfolding a graph into a tree */
  // igraph_unfold_tree

  /* Other Operations */
  // igraph_density
  // igraph_reciprocity
  // igraph_diversity
  // igraph_is_mutual
  // igraph_avg_nearest_neighbor_degree
  // igraph_get_adjacency
  // igraph_get_stochastic
  // Skipped igraph_get_stochastic_sparsemat
  // igraph_get_edgelist
  // igraph_contract_vertices

  /* Graph visitors */

  /* Breadth-first search
   * bfs(0, [](int vid, int pred, int succ, int rank, int dist) -> bool {
   *   // ...
   *   // Returns false to continue BFS. Details in doc of igraph_dfshandler_t
   *   return false;
   * });
   * */
  template <typename Function>
  void bfs(int root, Function callback = bfs_default_callback,
           Mode mode = Mode::Out, bool unreachable = true,
           const VectorView &restricted = VectorView::None(),
           Vector *order = nullptr, Vector *rank = nullptr,
           Vector *father = nullptr, Vector *pred = nullptr,
           Vector *succ = nullptr, Vector *dist = nullptr) {
    SafeCall(igraph_bfs(
        ptr(), root, NULL, static_cast<igraph_neimode_t>(mode), unreachable,
        restricted.ptr(), order ? order->ptr() : NULL,
        rank ? rank->ptr() : NULL, father ? father->ptr() : NULL,
        pred ? pred->ptr() : NULL, succ ? succ->ptr() : NULL,
        dist ? dist->ptr() : NULL,
        [](const igraph_t *, igraph_integer_t vid, igraph_integer_t pred,
           igraph_integer_t succ, igraph_integer_t rank, igraph_integer_t dist,
           void *extra) -> igraph_bool_t {
          auto cb = *reinterpret_cast<Function *>(extra);
          return static_cast<igraph_bool_t>(cb(vid, pred, succ, rank, dist));
        },
        reinterpret_cast<void *>(&callback)));
  }

  template <typename Function>
  void bfs(Vector roots, Function callback, Mode mode = Mode::Out,
           bool unreachable = true,
           const VectorView &restricted = VectorView::None()) {
    SafeCall(igraph_bfs(
        ptr(), NULL, roots.ptr(), static_cast<igraph_neimode_t>(mode),
        static_cast<igraph_bool_t>(unreachable), restricted.ptr(), NULL, NULL,
        NULL, NULL, NULL, NULL,
        [](const igraph_t *, igraph_integer_t vid, igraph_integer_t pred,
           igraph_integer_t succ, igraph_integer_t rank, igraph_integer_t dist,
           void *extra) -> igraph_bool_t {
          auto cb = *reinterpret_cast<Function *>(extra);
          return static_cast<igraph_bool_t>(cb(vid, pred, succ, rank, dist));
        },
        reinterpret_cast<void *>(&callback)));
  }

  /* Depth-first search
   * bfs(0, [](int vid, int dist) -> bool {
   *   // ...
   *   // Returns false to continue BFS. Details in doc of igraph_dfshandler_t
   *   return false;
   * });
   * */
  template <typename Function>
  void dfs(int root, Function in_callback, Mode mode = Mode::Out,
           bool unreachable = true, Vector *order = nullptr,
           Vector *order_out = nullptr, Vector *father = nullptr,
           Vector *dist = nullptr) {
    SafeCall(igraph_dfs(
        ptr(), root, static_cast<igraph_neimode_t>(mode), unreachable,
        order ? order->ptr() : NULL, order_out ? order_out->ptr() : NULL,
        father ? father->ptr() : NULL, dist ? dist->ptr() : NULL,
        [](const igraph_t *, igraph_integer_t vid, igraph_integer_t dist,
           void *extra) -> igraph_bool_t {
          auto cb = *reinterpret_cast<Function *>(extra);
          return static_cast<igraph_bool_t>(cb(vid, dist));
        },
        NULL, reinterpret_cast<void *>(&in_callback)));
  }

  template <typename Function, typename Function2>
  void dfs(int root, Function in_callback, Function2 out_callback,
           Mode mode = Mode::Out, bool unreachable = true,
           Vector *order = nullptr, Vector *order_out = nullptr,
           Vector *father = nullptr, Vector *dist = nullptr) {
    struct callbacks {
      callbacks(Function &icb, Function2 &ocb) : in_cb(icb), out_cb(ocb) {}
      Function &in_cb;
      Function2 &out_cb;
    } cb(in_callback, out_callback);
    SafeCall(igraph_dfs(
        ptr(), root, static_cast<igraph_neimode_t>(mode), unreachable,
        order ? order->ptr() : NULL, order_out ? order_out->ptr() : NULL,
        father ? father->ptr() : NULL, dist ? dist->ptr() : NULL,
        [](const igraph_t *, igraph_integer_t vid, igraph_integer_t dist,
           void *extra) -> igraph_bool_t {
          auto cb = *reinterpret_cast<callbacks *>(extra);
          return static_cast<igraph_bool_t>(cb.in_cb(vid, dist));
        },
        [](const igraph_t *, igraph_integer_t vid, igraph_integer_t dist,
           void *extra) -> igraph_bool_t {
          auto cb = *reinterpret_cast<callbacks *>(extra);
          return static_cast<igraph_bool_t>(cb.out_cb(vid, dist));
        },
        reinterpret_cast<void *>(&cb)));
  }

  /*  Cliques and Independent Vertex Sets */
  // Skipped igraph_cliques
  // Skipped igraph_largest_cliques
  // Skipped igraph_maximal_cliques
  // igraph_maximal_cliques_count
  // igraph_clique_number

  /* Independent Vertex Sets */
  // Skipped igraph_independent_vertex_sets
  // Skipped igraph_largest_independent_vertex_sets
  // Skipped igraph_maximal_independent_vertex_sets
  // igraph_independence_number

  /* Graph Isomorphism */
  // igraph_permute_vertices
  // igraph_isomorphic
  // igraph_subisomorphic
  // igraph_canonical_permutation
  // igraph_isomorphic_bliss
  // igraph_automorphisms
  // Skipped igraph_isomorphic_vf2
  // Skipped igraph_count_isomorphisms_vf2
  // Skipped igraph_get_isomorphisms_vf2
  // Skipped igraph_isomorphic_function_vf2
  // Skipped igraph_subisomorphic_vf2
  // Skipped igraph_count_subisomorphisms_vf2
  // Skipped igraph_get_subisomorphisms_vf2
  // Skipped igraph_subisomorphic_function_vf2
  // Skipped igraph_subisomorphic_lad
  // igraph_isomorphic_34
  // igraph_isoclass
  // igraph_isoclass_subgraph
  // igraph_isoclass_create

  /* Graph Motifs, Dyad Census and Triad Census */
  // igraph_dyad_census
  // igraph_triad_census
  // igraph_motifs_randesu
  // igraph_motifs_randesu_no
  // igraph_motifs_randesu_estimate
  // igraph_motifs_randesu_callback

  /* Generating Layouts for Graph Drawing */
  // igraph_layout_random
  // igraph_layout_circle
  // igraph_layout_star
  // igraph_layout_grid
  // igraph_layout_graphopt
  // Skipped igraph_layout_bipartite
  // Skipped igraph_layout_drl
  // Skipped igraph_layout_drl_3d
  // igraph_layout_fruchterman_reingold
  // igraph_layout_kamada_kawai
  // igraph_layout_mds
  // igraph_layout_grid_fruchterman_reingold
  // igraph_layout_lgl
  // igraph_layout_reingold_tilford
  // igraph_layout_reingold_tilford_circular
  // igraph_layout_sugiyama
  // igraph_layout_random_3d
  // igraph_layout_sphere
  // igraph_layout_grid_3d
  // igraph_layout_fruchterman_reingold_3d
  // igraph_layout_kamada_kawai_3d
  // Skipped igraph_layout_merge_dla

  /* Reading and Writing Graphs from and to Files */
  static Graph ReadEdgelist(FILE *instream, int n = 0,
                            bool directed = Undirected);
  static Graph ReadEdgelist(std::string filename, int n = 0,
                            bool directed = Undirected);
  // igraph_read_graph_edgelist
  // igraph_write_graph_edgelist
  // Skipped igraph_read_graph_ncol
  // igraph_write_graph_ncol
  // igraph_read_graph_lgl
  // igraph_write_graph_lgl(const
  // Skipped igraph_read_graph_dimacs
  // igraph_write_graph_dimacs
  // igraph_read_graph_graphdb
  // igraph_read_graph_graphml
  // igraph_write_graph_graphml
  // igraph_read_graph_gml
  // igraph_write_graph_gml
  static Graph ReadPajek(FILE *instream);
  static Graph ReadPajek(std::string filename);
  // igraph_read_graph_dl
  // igraph_write_graph_dot

  /* Maximum Flows, Minimum Cuts and related measures */
  // igraph_maxflow
  // igraph_maxflow_value
  // igraph_dominator_tree
  // igraph_st_mincut
  // igraph_st_mincut_value
  // Skipped igraph_all_st_cuts
  // Skipped igraph_all_st_mincuts
  // igraph_mincut
  // igraph_mincut_value

  /* Connectivity */
  // igraph_st_edge_connectivity(const
  // igraph_edge_connectivity
  // igraph_st_vertex_connectivity
  // igraph_vertex_connectivity

  /* Edge- and Vertex-Disjoint Paths */
  // igraph_edge_disjoint_paths
  // igraph_vertex_disjoint_paths

  /* Graph Adhesion and Cohesion */
  // igraph_adhesion
  // igraph_cohesion
  // Skipped igraph_cohesive_blocks

  /* Vertex separators */
  // igraph_is_separator
  // igraph_is_minimal_separator
  // Skipped igraph_all_minimal_st_separators
  // Skipped igraph_minimum_size_separators

  /* Detecting Community Structure */
  // igraph_modularity
  // igraph_community_optimal_modularity
  // igraph_community_to_membership
  // igraph_reindex_membership
  // igraph_compare_communities
  // igraph_split_join_distance
  // igraph_community_spinglass
  // igraph_community_spinglass_single
  // Skipped igraph_community_leading_eigenvector
  // igraph_le_community_to_membership
  // igraph_community_walktrap
  // igraph_community_edge_betweenness
  // igraph_community_eb_get_merges
  // igraph_community_fastgreedy
  // igraph_community_multilevel
  // igraph_community_label_propagation
  // igraph_community_infomap

  /* Graphlets */
  // Skipped igraph_graphlets
  // Skipped igraph_graphlets_candidate_basis
  // Skipped igraph_graphlets_project

  /* Hierarchical random graphs */
  // No hierarchical representation implemented

  /* Spectral Coarse Graining */
  // Skipped

  igraph_t *ptr() { return &graph_; }
  const igraph_t *ptr() const { return &graph_; }

  Graph(const igraph_t *graph) {
    graph_ = *graph;
    disown();
  }

 protected:
  Graph(const igraph_t &graph, bool owner = true);

 private:
  void disown() { owner_ = false; }
  bool owner() const { return owner_; }

  igraph_t graph_;
  bool owner_ = true;
};

template <>
struct TypeMapper<Graph> {
  typedef igraph_t type;

  static Graph Build(const type &graph, bool owner = false) {
    Graph g = Graph(graph);
    if (!owner) g.disown();
    return g;
  }
};

template <typename Iterator, typename>
Graph::Graph(Iterator edges_begin, Iterator edges_end, long int vertices,
             bool directed) {
  Vector vector(edges_begin, edges_end);
  SafeCall(igraph_create(ptr(), vector.ptr(), vertices, directed));
}

template <typename... Args, typename>
Graph Graph::LCF(int vertices, Args... args) {
  igraph_t graph;
  SafeCall(igraph_lcf(&graph, vertices, args..., 0));
  return Graph(graph);
}

}  // namespace igraph

#endif  // IGRAPHPP_GRAPH_HPP_
