#include <iostream>

#include <igraphpp/igraph.hpp>

int main(void) {
  using igraph::Graph;
  using std::cout;

  igraph::SetSeed(17);
  Graph g = Graph::ErdosRenyiGame(50, 0.1, igraph::Undirected, igraph::NoLoops);
  cout << "Vertices : " << g.vcount() << "\n";
  cout << "Edges    : " << g.ecount() << "\n";
  cout << "Directed : " << (g.is_directed()? "yes": "no") << "\n";
  cout << "Connected: " << (g.is_connected()? "yes": "no") << "\n";
  cout << "Diameter : " << g.diameter() << "\n";
  return 0;
}
