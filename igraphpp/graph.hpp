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
  Graph &add_edges(const Vector &edges);
  Graph &add_edges(std::initializer_list<double> edges);
  Graph &add_vertices(int number_of_vertices);
  void delete_edges(const EdgeSelector &edges);
  void delete_edges(std::initializer_list<double> edges);
  void delete_vertices(const VertexSelector &vertices);
  void delete_vertices(std::initializer_list<double> vertices);

  bool is_connected(Connectedness mode = WeaklyConnected) const;

  int diameter(Directedness directed = Directed, bool unconnected = true) const;

  double AveragePathLength(Directedness directed = Directed,
                           bool unconnected = true) const;

  /* Deterministic graph generators */
  static Graph AdjacencyMatrix(Matrix &adjmatrix,
                               AdjacencyMatrixMode mode = AdjacencyDirected);
  // Skipped igraph_weighted_adjacency
  // Skipped  igraph_adjlist
  static Graph Star(int vertices, StarMode mode = StarOut,
                    int center_vertex = 0);
  static Graph Lattice(const Vector &dimension, int nei = 1,
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
  static Graph Barabasi(int vertices, double power, const Vector &outseq,
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
  static Graph StaticFitness(int vertices, VectorView &fitness,
                             Loops loops = NoLoops, bool multiple = false);
  static Graph StaticFitness(int vertices, VectorView &fitness_out,
                             VectorView &fitness_in, Loops loops = NoLoops,
                             bool multiple = false);
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
