
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("EdgeIterator", "[EdgeIterator]") {
  using igraph::EdgeIterator;
  using igraph::EdgeSelector;
  using igraph::Graph;

  Graph g(10);
  g.AddEdge(0, 1).AddEdge(0, 2).AddEdge(0, 3).AddEdge(1, 2).AddEdge(4, 6);

  EdgeIterator it(g, EdgeSelector::All());
  CHECK(it.size() == 5);

  long int x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }
  CHECK(x == 5);

  x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }
  CHECK(x == 5);
}
