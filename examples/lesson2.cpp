#include <iostream>

#include <igraphpp/igraph.hpp>

int main(void) {
  using igraph::Graph;
  using igraph::Vector;
  using std::cout;

  Vector dim{{30, 30}};
  Graph graph = Graph::Lattice(dim, 0, igraph::Undirected, igraph::NotMutual,
                               igraph::Circular);

  igraph::SetSeed(42);

  Vector edges(20);
  for (long int i = 0; i < edges.size(); i++) {
    edges[i] = rand() % graph.vcount();
  }

  double avg_path = graph.average_path_length(igraph::Undirected);
  cout << "Average path length (lattice):            " << avg_path << "\n";

  graph.add_edges(edges);
  avg_path = graph.average_path_length(igraph::Undirected);
  cout << "Average path length (randomized lattice): " << avg_path << "\n";

  return 0;
}
