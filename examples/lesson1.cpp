#include <iostream>

#include <igraphpp/igraph.hpp>

int main(void) {
  using igraph::Graph;
  using std::cout;

  igraph::SetSeed(42);
  Graph g = Graph::ErdosRenyiGame(1000, 5.0 / 1000);

  cout << "Vertices : " << g.vcount() << "\n";
  cout << "Edges    : " << g.ecount() << "\n";
  cout << "Directed : " << (g.is_directed() ? "yes" : "no") << "\n";
  cout << "Connected: " << (g.is_connected() ? "yes" : "no") << "\n";
  cout << "Diameter : " << g.diameter() << "\n";
  return 0;
}
