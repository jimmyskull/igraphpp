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
  Graph &add_edge(int from, int to);
  Graph &add_edges(const VectorView &edges);
  Graph &add_edges(std::initializer_list<double> edges);
  Graph &add_vertices(int number_of_vertices);
  void delete_edges(const EdgeSelector &edges);
  void delete_edges(std::initializer_list<double> edges);
  void delete_vertices(const VertexSelector &vertices);
  void delete_vertices(std::initializer_list<double> vertices);

  /* Deterministic graph generators */
  static Graph AdjacencyMatrix(Matrix &adjmatrix,
                               AdjacencyMatrixMode mode = AdjacencyDirected);
  // Skipped igraph_weighted_adjacency (attributes has not been implemented)
  // Skipped  igraph_adjlist (igraph_adjlist_t has not been implemented)
  static Graph Star(int vertices, StarMode mode = StarOut,
                    int center_vertex = 0);
  static Graph Lattice(const VectorView &dimension, int nei = 1,
                       Directedness dir = Undirected,
                       Mutuality mutual = NotMutual,
                       Periodicity periodicity = NotPeriodic);
  static Graph Ring(int vertices, Directedness dir = Undirected,
                    Mutuality mutual = NotMutual, bool circular = true);
  static Graph Tree(int vertices, int children = 2, TreeMode mode = TreeOut);
  static Graph Full(int vertices, Directedness dir = Undirected,
                    Loops loops = NoLoops);
  static Graph FullCitation(int vertices, Directedness dir = Directed);
  static Graph Famous(const char *name);
  template <typename... Args, typename = std::enable_if_t<util::all_args(
                                  std::is_same<Args, int>::value...)>>
  static Graph LCF(int vertices, Args... args);
  static Graph LCF(int vertices, const VectorView &shifts, int repeats);
  static Graph LCF(int vertices, std::initializer_list<double> shifts,
                   int repeats);
  static Graph Atlas(int number);
  static Graph deBruijn(int m, int n);
  static Graph Kautz(int m, int n);
  static Graph ExtendedChordalRing(int vertices, const Matrix &W);
  void connect_neighborhood(int order, NeighborMode mode = Out);

  /* Randomized graph generators */
  static Graph GRG(int vertices, double radius, bool torus);
  static Graph Barabasi(int vertices, double power = 1, int m = 2,
                        bool outpref = false, double A = 1,
                        Directedness dir = Directed,
                        BarabasiAlgorithm algo = BarabasiPSumTree,
                        const Graph &start_graph = Graph());
  static Graph Barabasi(int vertices, double power, const VectorView &outseq,
                        bool outpref = false, double A = 1,
                        Directedness dir = Directed,
                        BarabasiAlgorithm algo = BarabasiPSumTree,
                        const Graph &start_graph = Graph());
  static Graph ErdosRenyi(int vertices, double prob,
                          Directedness dir = Undirected, Loops loops = NoLoops);
  static Graph ErdosRenyi(int vertices, int edges,
                          Directedness dir = Undirected, Loops loops = NoLoops);
  static Graph WattsStrogatz(int dim, int size, int nei, double p,
                             Loops loops = NoLoops, bool multiple = false);
  // Rewire uniformly
  void rewire_edges(double prob, Loops loops = NoLoops, bool multiple = false);
  static Graph
  DegreeSequence(const VectorView &degrees,
                 DegreeSequenceMethod method = DegreeSequenceSimple);
  static Graph
  DegreeSequence(const VectorView &out_deq, const VectorView &in_deq,
                 DegreeSequenceMethod method = DegreeSequenceSimple);
  static Graph kRegular(int vertices, int k, Directedness dir = Directed,
                        bool multiple = false);
  static Graph StaticFitness(int vertices, const VectorView &fitness,
                             Loops loops = NoLoops, bool multiple = false);
  static Graph StaticFitness(int vertices, const VectorView &fitness_out,
                             const VectorView &fitness_in,
                             Loops loops = NoLoops, bool multiple = false);
  static Graph StaticPowerLaw(int vertices, int edges, double exponent_out,
                              double exponent_in = -1.0, Loops loops = NoLoops,
                              bool multiple = false,
                              bool finite_size_correlation = true);
  static Graph ForestFire(int vertices, double fw_prob, double bw_factor = 1.0,
                          int pambs = 1, Directedness dir = Undirected);
  // Rewire keeping the degree distribution
  void rewire(int trials, RewiringMode mode = RewiringSimple);
  static Graph GrowingRandom(int vertices, int m = 1,
                             Directedness dir = Directed,
                             bool citation = false);
  static Graph CallawayTraits(int vertices, int types, int edges_per_step,
                              const VectorView &type_dist,
                              const Matrix &pref_matrix,
                              Directedness dir = Undirected);
  static Graph Establishment(int vertices, int types, int trials_per_step,
                             const VectorView &type_dist,
                             const Matrix &pref_matrix,
                             Directedness dir = Undirected);
  static Graph Preference(int vertices, int types, const VectorView &type_dist,
                          bool fixed_sizes, const Matrix &pref_matrix,
                          Directedness dir = Undirected, Loops loops = NoLoops);
  static Graph AsymmetricPreference(int vertices, int types,
                                    const Matrix &type_dist_matrix,
                                    const Matrix &pref_matrix,
                                    Loops loops = NoLoops);
  static Graph RecentDegree(int vertices, double power, int window,
                            bool outpref = false, double zero_appeal = 1.0,
                            Directedness dir = Undirected);
  static Graph RecentDegree(int vertices, double power, int window, int m,
                            const VectorView &outseq, bool outpref = false,
                            double zero_appeal = 1.0,
                            Directedness dir = Undirected);
  static Graph BarabasiAging(int vertices, int m, bool outpref, double pa_exp,
                             double aging_exp, int aging_bin,
                             double zero_deg_appeal, double zero_age_appeal,
                             double deg_coef, double age_coef,
                             Directedness dir = Undirected);
  static Graph BarabasiAging(int vertices, const VectorView &outseq,
                             bool outpref, double pa_exp, double aging_exp,
                             int aging_bin, double zero_deg_appeal,
                             double zero_age_appeal, double deg_coef,
                             double age_coef, Directedness dir = Undirected);
  static Graph RecentDegreeAging(int vertices, int m, bool outpref,
                                 double pa_exp, double aging_exp, int aging_bin,
                                 double time_window, double zero_appeal,
                                 Directedness dir = Undirected);
  static Graph RecentDegreeAging(int vertices, const VectorView &outseq,
                                 bool outpref, double pa_exp, double aging_exp,
                                 int aging_bin, double time_window,
                                 double zero_appeal,
                                 Directedness dir = Undirected);
  static Graph CitedType(int vertices, const VectorView &types,
                         const VectorView &pref, int edges_per_step,
                         Directedness dir = Undirected);
  // Skipped igraph_sbm_game (igraph_vector_int_t has not been implemented)

  /* Graph, vertex, and edge attributes */

  /* Structural properties of graphs */
  /* Basic properties */
  bool are_connected(int v1, int v2) const;

  /* Shortest path related functions */
  Matrix shortest_paths(const VertexSelector &from, const VertexSelector &to,
                        NeighborMode mode = Out) const;
  double shortest_paths(int from, int to, NeighborMode mode = Out) const;
  Matrix shortest_paths_dijkstra(const VertexSelector &from,
                                 const VertexSelector &to,
                                 const VectorView &weights,
                                 NeighborMode mode = Out) const;
  double shortest_paths_dijkstra(int from, int to, const VectorView &weights,
                                 NeighborMode mode = Out) const;
  Matrix shortest_paths_bellman_ford(const VertexSelector &from,
                                     const VertexSelector &to,
                                     const VectorView &weights,
                                     NeighborMode mode = Out) const;
  double shortest_paths_bellman_ford(int from, int to,
                                     const VectorView &weights,
                                     NeighborMode mode = Out) const;
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
  double average_path_length(Directedness dir = Directed,
                             bool unconnected = true) const;
  Vector path_length_hist(Directedness dir = Directed,
                          double *unconnected = NULL) const;
  int diameter(int *source_vertex = NULL, int *target_vertex = NULL,
               Directedness dir = Directed, bool unconnected = true) const;
  int diameter(VectorView &path, Directedness dir = Directed,
               bool unconnected = true) const;
  double diameter_dijkstra(const VectorView &weights, int *source_vertex = NULL,
                           int *target_vertex = NULL,
                           Directedness dir = Directed,
                           bool unconnected = true) const;
  double diameter_dijkstra(const VectorView &weights, VectorView &path,
                           Directedness dir = Directed,
                           bool unconnected = true) const;
  int girth() const;
  int girth(Vector &circle) const;
  Vector eccentricity(const VertexSelector &vids,
                      NeighborMode mode = Out) const;
  Vector eccentricity(int vertex, NeighborMode mode = Out) const;
  double radius(NeighborMode mode = Out) const;

  bool is_connected(Connectedness mode = WeaklyConnected) const;

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
