
#include <vector>

#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("EdgeSelector", "[EdgeSelector]") {
  using igraph::EdgeSelector;
  using igraph::Vector;
  using igraph::Graph;

  Graph g(10);
  g.AddEdge(0, 1).AddEdge(0, 2).AddEdge(2, 3);

  EdgeSelector vall = EdgeSelector::All();
  CHECK(vall.size(g) == 3);
  CHECK(vall.is_all());

  EdgeSelector vinc = EdgeSelector::Incident(0);
  CHECK(vinc.size(g) == 2);
  CHECK_FALSE(vinc.is_all());

  EdgeSelector vnone = EdgeSelector::None();
  CHECK(vnone.size(g) == 0);
  CHECK_FALSE(vnone.is_all());

  EdgeSelector vsingle = EdgeSelector::Single(0);
  CHECK(vsingle.size(g) == 1);
  CHECK_FALSE(vsingle.is_all());

  Vector v(1, 10);
  CHECK(v.size() == 10);
  EdgeSelector vvector = EdgeSelector::FromVector(v);
  EdgeSelector vseq = EdgeSelector::Sequence(0, 3);
  EdgeSelector copy(vseq);
  copy = copy;
  CHECK_FALSE(copy.is_all());
  CHECK(copy.size(g) == 3);

  EdgeSelector pairs = EdgeSelector::Pairs(igraph::Undirected, 0, 2);
  CHECK(pairs.size(g) == 1);
  CHECK_FALSE(pairs.is_all());

  EdgeSelector vseqall = EdgeSelector::Sequence(0, 9);
  CHECK_FALSE(vseqall.is_all());

  EdgeSelector einit({1, 2, 2, 1});
  CHECK(einit.size(g) == 4);

  std::vector<double> vd = {0, 1, 1, 2, 3, 4};
  EdgeSelector einit2(vd.begin(), vd.end());
  CHECK(einit2.size(g) == 6);

  EdgeSelector einit3({6, 7, 8, 9, 0, 1, 0, 2});
  CHECK(einit3.size(g) == 8);
}
