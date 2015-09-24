
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("VertexIterator", "[VertexIterator]") {
  using igraph::VertexIterator;
  using igraph::VertexSelector;
  // using igraph::Vector;
  using igraph::Graph;

  Graph g(10);

  VertexIterator it(g, VertexSelector::All());

  long int x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }

  x = 0;
  for (auto i : it) {
    CHECK(x == i);
    ++x;
  }
}
