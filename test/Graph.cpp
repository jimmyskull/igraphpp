
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

#define N 10

TEST_CASE("Graph — contructors", "[Graph]") {
  using igraph::Graph;
  using igraph::VectorView;

  Graph empty;
  CHECK(empty.vcount() == 0);
  CHECK(empty.ecount() == 0);

  Graph empty_undir(N, igraph::Undirected);
  CHECK(empty_undir.vcount() == N);
  CHECK(empty_undir.ecount() == 0);
  CHECK_FALSE(empty_undir.is_directed());

  Graph empty_dir(N, igraph::Directed);
  CHECK(empty_dir.vcount() == N);
  CHECK(empty_dir.ecount() == 0);
  CHECK(empty_dir.is_directed());

  Graph initlist({1, 3, 1, 6});
  CHECK(initlist.vcount() == 7);
  CHECK(initlist.ecount() == 2);

  std::vector<double> vd = {1, 3, 1, 6};
  Graph initvec(vd.begin(), vd.end());
  CHECK(initvec.vcount() == 7);
  CHECK(initvec.ecount() == 2);

  Graph gmv(N);
  CHECK_FALSE(gmv.is_directed());
  Graph cpy(gmv);
  Graph gmv1(std::move(gmv));
  Graph gmv2(std::move(gmv1));

  Graph assign(N);
  Graph assign2 = assign;
  Graph assign3 = std::move(assign);
}

TEST_CASE("Graph — basic query operations", "[Graph]") {
  using igraph::Graph;
  using igraph::Edge;
  using igraph::Vector;
  using igraph::VertexSelector;

  Graph g({0, 1, 0, 2, 1, 4, 1, 2, 3, 1});
  CHECK(g.vcount() == 5);
  CHECK(g.ecount() == 5);
  Edge e = g.edge(0);
  CHECK(e.first == 0);
  CHECK(e.second == 1);
  e = g.edge(2);
  CHECK(e.first == 1);
  CHECK(e.second == 4);
  CHECK(g.eid(0, 2) == 1);
  CHECK(g.eid(0, 4) == -1);

  Vector pairs{{0, 2, 2, 1}};
  Vector eids = g.pairs_eids(pairs);
  CHECK(eids.size() == 2);
  CHECK(eids[0] == 1);
  CHECK(eids[1] == 3);

  eids = g.pairs_eids({0, 2, 2, 1});
  CHECK(eids.size() == 2);
  CHECK(eids[0] == 1);
  CHECK(eids[1] == 3);

  Vector path{{0, 1, 4, 3, 1}};
  Vector path_eids = g.path_eids(path);
  CHECK(path_eids.size() == 4);
  Vector check_path{{0, 2, -1, 4}};
  CHECK(path_eids == check_path);

  path_eids = g.path_eids({0, 1, 4, 3, 1});
  CHECK(path_eids.size() == 4);
  CHECK(path_eids == check_path);

  eids = g.eids(pairs, path);
  CHECK(eids.size() == 6);
  check_path = Vector{{1, 3, 0, 2, -1, 4}};
  CHECK(eids == check_path);

  eids = g.eids({0, 2, 2, 1}, {0, 1, 4, 3, 1});
  CHECK(eids == check_path);

  CHECK(g.neighbors(0) == Vector({1, 2}));

  Vector ic = g.incident(2);
  CHECK(g.incident(2) == Vector({1, 3}));

  VertexSelector sel{{0, 1, 3}};
  CHECK(g.degree(sel) == Vector({2, 4, 1}));

  CHECK(g.degree({0, 1, 3, 2}) == Vector({2, 4, 1, 2}));
}
