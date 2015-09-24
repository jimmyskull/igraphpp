
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("VertexIterator", "[VertexIterator]") {
  using igraph::VertexIterator;
  using igraph::VertexSelector;
  using igraph::Graph;

  Graph g(10);

  VertexIterator it(g, VertexSelector::All());
  CHECK(it.size() == 10);

  long int x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }
  CHECK(x == 10);

  x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }
  CHECK(x == 10);
}
