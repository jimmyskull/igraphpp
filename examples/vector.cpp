#include <iostream>

#include <igraphpp/igraph.hpp>

int main(void) {
  using igraph::Graph;
  using igraph::VectorView;
  using std::cout;

  double edges[10] = {0, 1, 0, 2, 1, 2, 4, 5, 6, 3};
  int N = 9;

  Graph g = Graph(N, VectorView(edges, sizeof(edges) / sizeof(0[edges])));

  cout << "Vertices : " << g.vcount() << "\n";
  cout << "Edges    : " << g.ecount() << "\n";
  return 0;
}
