#include <iostream>

#include <igraphpp/igraph.hpp>

#include <time.h>

int main(void) {
  using igraph::Graph;
  using std::cout;

  igraph::SetSeed(42);
  float startTime = (float)clock() / CLOCKS_PER_SEC;
  Graph g = Graph::ErdosRenyiGame(1000, 5.0 / 1000);
  float endTime = (float)clock() / CLOCKS_PER_SEC;
  float timeElapsed = endTime - startTime;
  cout << "Time: " << timeElapsed << "\n";

  cout << "Vertices : " << g.vcount() << "\n";
  cout << "Edges    : " << g.ecount() << "\n";
  cout << "Directed : " << (g.is_directed() ? "yes" : "no") << "\n";
  cout << "Connected: " << (g.is_connected() ? "yes" : "no") << "\n";
  cout << "Diameter : " << g.diameter() << "\n";
  return 0;
}
