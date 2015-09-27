
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

TEST_CASE("Graph — basic properties", "[Graph]") {
  using igraph::Graph;

  Graph g{{0, 1, 0, 2}};
  CHECK(g.are_connected(0, 1));
  CHECK_FALSE(g.are_connected(1, 2));
}

TEST_CASE("Graph — shortest path related functions", "[Graph]") {
  using igraph::Graph;
  using igraph::Matrix;
  using igraph::VertexSelector;
  using igraph::Vector;

  Graph g{{0, 1, 0, 2, 1, 3}};
  double shortestpath = g.shortest_paths(0, 3);
  CHECK(shortestpath == 2);

  Matrix spmat =
      g.shortest_paths(VertexSelector::Single(0), VertexSelector::Single(3));
  CHECK(spmat(0, 0) == 2);

  Matrix mat = g.shortest_paths(VertexSelector::All(), VertexSelector::All());
  CHECK(mat.nrow() == 4);
  CHECK(mat.ncol() == 4);
  CHECK(mat == Matrix({0, 1, 1, 2, 1, 0, 2, 1, 1, 2, 0, 3, 2, 1, 3, 0}, 4));

  g.add_edge(0, 3);
  CHECK(g.shortest_paths_dijkstra(0, 3, Vector({1, 1, 1, 5})) == 2);
  CHECK(g.shortest_paths_dijkstra(0, 3, Vector({1, 1, 1, 1})) == 1);
  CHECK(g.shortest_paths_bellman_ford(0, 3, Vector({1, 1, 1, 5})) == 2);
  CHECK(g.shortest_paths_bellman_ford(0, 3, Vector({2, 1, 1, 0})) == 0);
  CHECK(g.shortest_paths_johnson(0, 3, Vector({1, 1, 1, 5})) == 2);
  CHECK(g.shortest_paths_johnson(0, 3, Vector({2, 1, 1, 0})) == 0);

  // Skipped igraph_get_shortest_paths
  // Skipped igraph_get_shortest_path
  // Skipped igraph_get_shortest_paths_dijkstra
  // Skipped igraph_get_shortest_path_dijkstra
  // Skipped igraph_get_all_shortest_paths
  // Skipped igraph_get_all_shortest_paths_dijkstra

  CHECK(g.average_path_length() == Approx(1.33333));
  Vector hist = g.path_length_hist();
  CHECK(hist == Vector({4, 2}));
  int pfrom, pto;
  CHECK(g.diameter(&pfrom, &pto) == 2);
  CHECK(pfrom == 1);
  CHECK(pto == 2);
  Vector spath;
  CHECK(g.diameter(spath) == 2);
  CHECK(spath == Vector({1, 0, 2}));
  CHECK(spath.head() == 1);
  CHECK(spath.tail() == 2);

  CHECK(g.diameter_dijkstra(Vector({0.5, 0.5, 1, 1})) == 1.5);

  Vector circle;
  CHECK(g.girth(circle) == 3);
  CHECK(circle.size() == 3);
  CHECK(circle == Vector({1, 0, 3}));

  CHECK(g.eccentricity(VertexSelector::All()) == Vector({1, 2, 2, 2}));
  CHECK(g.eccentricity(1) == 2);
  CHECK(g.radius() == 1);
}

TEST_CASE("Graph — neighborhood of a vertex", "[Graph]") {
  using igraph::Graph;
  using igraph::Vector;

  Graph g{{0, 1, 0, 2, 0, 3, 0, 4, 2, 5, 2, 6, 3, 5}};
  REQUIRE(g.ecount() == 7);
  REQUIRE(g.vcount() == 7);

  CHECK(g.neighborhood_size(0) == 5);
  CHECK(g.neighborhood_size({{1, 5}}) == Vector({2, 3}));
}

TEST_CASE("Graph — graph components", "[Graph]") {
  using igraph::Graph;
  using igraph::Vector;

  Graph g{{0, 1, 0, 2, 0, 3, 0, 4, 2, 5, 2, 6, 8, 9, 8, 10}};
  REQUIRE(g.ecount() == 8);
  REQUIRE(g.vcount() == 11);

  CHECK(g.subcomponent(0).size() == 7);
  CHECK(g.subcomponent(10).size() == 3);
  Graph sub(g.induced_subgraph(g.subcomponent(0)));
  CHECK(sub.vcount() == 7);
  CHECK(g.subgraph_edges({{1, 4}}).diameter() == 2);

  CHECK(g.clusters().number_of_clusters == 3);

  CHECK_FALSE(g.is_connected());

  CHECK(g.articulation_points().size() == 3);
}

TEST_CASE("Graph — degree sequences", "[Graph]") {
  using igraph::Graph;

  CHECK_FALSE(Graph::is_degree_sequence({{0, 2, 3, 0, 4, 3, 1, 3, 4, 2}},
                                        {{0, 3, 1, 2, 2, 4, 4, 1, 3, 1}}));

  CHECK(Graph::is_graphical_degree_sequence({{1, 3, 2, 1, 3, 4, 3, 3, 1, 3}},
                                            {{4, 1, 2, 3, 2, 3, 2, 3, 2, 2}}));
}

TEST_CASE("Graph — centrality measures", "[Graph]") {
  using igraph::Graph;
  using igraph::VertexSelector;
  using igraph::Vector;
  using igraph::PageRank;

  Graph g = Graph::Star(10);
  Vector closeness = g.closeness(VertexSelector::All());
  CHECK(closeness[0] == Approx(0.111111));
  for (int i = 1; i < 10; ++i)
    CHECK(closeness[i] == Approx(0.0111111));
  CHECK(g.closeness(0, igraph::Out, Vector(1, g.ecount())) ==
        Approx(0.0222222222));

  g.add_edges({1, 2, 2, 3, 4, 7});
  Vector betweenness = g.betweenness(VertexSelector::All());
  CHECK(betweenness == Vector({0, 0, 1, 0, 0, 0, 0, 0, 0, 0}));
  CHECK(g.betweenness(2) == 1);

  betweenness = g.edge_betweenness();
  CHECK(betweenness == Vector({1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1}));

  PageRank pr = g.pagerank();
  CHECK(pr.scores[4] == Approx(0.0758368));
}
