
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("VectorPtr<Graph>", "[VectorPtr]") {
  using igraph::VectorPtr;
  using igraph::Graph;

  Graph g1(1), g2(2), g3(3);

  VectorPtr<Graph> generic;
  CHECK(generic.size() == 0);
  generic.push_back(g1);
  generic.push_back(g2);
  CHECK(generic.size() == 2);
  CHECK(generic[0].vcount() == g1.vcount());
  generic.clear();
  CHECK(generic.size() == 0);

  VectorPtr<Graph> vector;
  for (long int i = 0; i < 10; ++i) {
    Graph g(i);
    vector.push_back(g);
    CHECK(Graph(vector[i]).vcount() == g.vcount());
  }
  generic.clear();
}

TEST_CASE("VectorPtr<VectorView>", "[VectorPtr]") {
  using igraph::VectorPtr;
  using igraph::VectorView;
  using igraph::Vector;

  const long int size = 5;
  double data[size] = {1, 2, 3, 4, 5};

  VectorView view(data, size);
  Vector seq(0.0, 20.0);

  VectorPtr<VectorView> vector;
  vector.push_back(view);
  vector.push_back(seq);

  CHECK(vector.size() == 2);
  CHECK(vector[0].size() == size);
  CHECK(vector[1].size() == 21);
}

