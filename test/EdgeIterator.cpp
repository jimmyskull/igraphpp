
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("EdgeIterator", "[EdgeIterator]") {
  using igraph::EdgeIterator;
  using igraph::EdgeSelector;
  using igraph::Graph;

  Graph g(10);
  g.add_edge(0, 1).add_edge(0, 2).add_edge(0, 3).add_edge(1, 2).add_edge(4, 6);

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
