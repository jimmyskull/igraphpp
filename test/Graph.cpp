
#include <cmath>

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

TEST_CASE("Graph — adding and deleting vertices and edges", "[Graph]") {
  using igraph::Graph;
  using igraph::Vector;
  using igraph::EdgeSelector;
  using igraph::VertexSelector;

  Graph g(10);
  g.add_edge(0, 1);
  CHECK(g.ecount() == 1);

  Vector elist{{0, 2, 0, 3, 0, 4}};
  g.add_edges(elist);
  CHECK(g.ecount() == 4);

  g.add_edges({1, 2, 2, 3, 2, 4, 4, 3});
  CHECK(g.ecount() == 8);

  g.add_vertices(10);
  CHECK(g.vcount() == 20);

  CHECK(g.eid(0, 2) == 1);
  g.delete_edges(EdgeSelector::Pairs({0, 2}));
  CHECK(g.eid(0, 2) == -1);
  CHECK(g.eid(0, 3) == 1);
  CHECK(g.eid(0, 4) == 2);
  g.delete_edges({1, 2});
  CHECK(g.eid(0, 3) == -1);
  CHECK(g.eid(0, 4) == -1);
  CHECK(g.eid(1, 2) == 1);

  CHECK(g.eid(1, 2) == 1);
  // All vertices higher than 1 are decreased by one
  g.delete_vertices(VertexSelector::Single(1));

  CHECK(g.pairs_eids({2, 3, 2, 1}) > Vector({-1, -1}));
  g.delete_vertices({2, 3, 2, 1});
  CHECK(g.pairs_eids({2, 3, 2, 1}) == Vector({-1, -1}));
}

TEST_CASE("Graph — deterministic generators", "[Graph]") {
  using igraph::Matrix;
  using igraph::Graph;
  using igraph::Vector;

  Matrix mat(N, N);
  mat(0, 1) = 1;
  mat(1, 0) = 1;
  mat(1, 3) = 1;
  mat(1, 2) = 1;
  Graph g = Graph::AdjacencyMatrix(mat, igraph::AdjacencyUndirected);
  CHECK(g.vcount() == 10);
  CHECK(g.ecount() == 3);

  Graph star = Graph::Star(10, igraph::StarOut, 3);
  CHECK(star.degree({3})[0] == 9);

  Graph lat = Graph::Lattice(Vector({3, 3}), 1);
  CHECK(lat.vcount() == 9);
  CHECK(lat.ecount() == 12);

  Graph ring = Graph::Ring(10);
  CHECK(ring.vcount() == 10);
  CHECK(ring.ecount() == 10);

  Graph tree = Graph::Tree(10);
  CHECK(tree.vcount() == 10);
  CHECK(tree.ecount() == 9);

  Graph full = Graph::Full(10);
  CHECK(full.vcount() == 10);
  CHECK(full.ecount() == 10 * 9 / 2);

  Graph full_citation = Graph::FullCitation(10);
  CHECK(full_citation.is_directed());
  CHECK(full_citation.vcount() == 10);
  CHECK(full_citation.ecount() == 10 * 9 / 2);

  Graph karate = Graph::Famous("Zachary");
  CHECK(karate.vcount() == 34);

  Graph lcf = Graph::LCF(12, 5, -5, 6);
  Vector shifts{{5, -5}};
  lcf = Graph::LCF(12, shifts, 6);
  lcf = Graph::LCF(12, {5, -5}, 6);

  Graph atlas = Graph::Atlas(1);
  Graph bruijn = Graph::deBruijn(10, 2);
  Graph kautz = Graph::Kautz(10, 2);

  Matrix W(2, 3);
  W.set_row(0, Vector{{3, 4, 8}});
  W.set_row(1, Vector{{12, 7, 11}});
  Graph ecr = Graph::ExtendedChordalRing(15, W);

  Graph gnei{{0, 1, 1, 2, 2, 3}};
  CHECK(gnei.ecount() == 3);
  gnei.connect_neighborhood(2);
  CHECK(gnei.ecount() == 5);
}

TEST_CASE("Graph — randomized graph generators", "[Graph]") {
  using igraph::Matrix;
  using igraph::Graph;
  using igraph::Vector;

  Graph grg = Graph::GRG(10, std::sqrt(2.0) / 2.0, true);
  CHECK(grg.vcount() == 10);
  CHECK(grg.ecount() == 10 * 9 / 2);

  Graph barabasi = Graph::Barabasi(100);

  Graph er = Graph::ErdosRenyi(100, 0.1);
  CHECK(er.vcount() == 100);

  Graph ws = Graph::WattsStrogatz(1, 100, 5, 0.01);
  CHECK(ws.vcount() == 100);
  ws.rewire_edges(0.4);

  Graph ds = Graph::DegreeSequence(Vector({1, 2, 3, 1, 1}));
  CHECK_FALSE(ds.is_directed());
  ds = Graph::DegreeSequence(Vector({1, 2, 3, 1, 1}), Vector({1, 2, 3, 1, 1}));
  CHECK(ds.is_directed());

  Graph reg = Graph::kRegular(10, 2);
  CHECK(reg.vcount() == 10);
  CHECK(reg.ecount() == 20);

  Graph spl = Graph::StaticPowerLaw(10, 20, 2.0);

  Graph ff = Graph::ForestFire(200, 0.37, 0.32 / 0.37);
  ff.rewire(100);

  Graph ct =
      Graph::CallawayTraits(200, 2, 1, Vector({1, 1}), Matrix({1, 0, 0, 1}, 2));
  Graph st =
      Graph::Establishment(200, 2, 2, Vector({1, 1}), Matrix({1, 0, 0, 1}, 2));

  Graph pref =
      Graph::Preference(20, 2, Vector({1, 1}), false, Matrix({1, 0, 0, 1}, 2));
  Graph apref = Graph::AsymmetricPreference(20, 2, Matrix({1, 1, 1, 1}, 2),
                                            Matrix({0, 1, 0, 0}, 2));

  // RecentDegree
  // BarabasiAging
  // RecentDegreeAging
  // CitedType
}
