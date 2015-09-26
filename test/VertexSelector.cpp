
#include <vector>

#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("VertexSelector", "[VertexSelector]") {
  using igraph::VertexSelector;
  using igraph::VectorView;
  using igraph::Vector;
  using igraph::Graph;

  Graph g(10);
  g.AddEdge(0, 1);

  VertexSelector vall = VertexSelector::All();
  CHECK(vall.size(g) == 10);
  CHECK(vall.is_all());

  VertexSelector vadj = VertexSelector::Adjacent(0);
  CHECK(vadj.size(g) == 1);
  CHECK_FALSE(vadj.is_all());

  VertexSelector vnadj = VertexSelector::NonAdjacent(0);
  CHECK(vnadj.size(g) == 9);
  CHECK_FALSE(vnadj.is_all());

  VertexSelector vnone = VertexSelector::None();
  CHECK(vnone.size(g) == 0);
  CHECK_FALSE(vnone.is_all());

  VertexSelector vsingle = VertexSelector::Single(0);
  CHECK(vsingle.size(g) == 1);
  CHECK_FALSE(vsingle.is_all());

  Vector v(1, 10);
  CHECK(v.size() == 10);
  VertexSelector vvector = VertexSelector::FromVector(v);
  VertexSelector vseq = VertexSelector::Sequence(1, 7);
  VertexSelector copy(vseq);
  copy = copy;
  CHECK_FALSE(copy.is_all());
  CHECK(copy.size(g) == 7);

  VertexSelector small = VertexSelector::Small(1, 2, 3);
  CHECK(small.size(g) == 3);
  CHECK_FALSE(small.is_all());

  VertexSelector vseqall = VertexSelector::Sequence(0, 9);
  CHECK_FALSE(vseqall.is_all());

  VertexSelector vinit({1, 2, 3});
  CHECK(vinit.size(g) == 3);

  std::vector<double> vd = {0, 1, 2, 3, 4};
  VertexSelector vinit2(vd.begin(), vd.end());
  CHECK(vinit2.size(g) == 5);

  VertexSelector vinit3({6, 7, 8, 9});
  CHECK(vinit3.size(g) == 4);
}
