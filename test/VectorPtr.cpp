
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

TEST_CASE("VectorPtr", "[VectorPtr]") {
  using igraph::VectorPtr;
  using igraph::Graph;

  VectorPtr<void *> generic;
  CHECK(generic.size() == 0);
  generic.push_back(NULL);
  generic.push_back(NULL);
  CHECK(generic.size() == 2);
  CHECK(generic[0] == NULL);
  generic.clear();
  CHECK(generic.size() == 0);

  using VectorGraph = VectorPtr<igraph_t *>;

  VectorGraph vector;
  for (long int i = 0; i < 10; ++i) {
    Graph g(i);
    vector.push_back(g.ptr());
    CHECK(Graph(vector[i]).vcount() == i);
  }
  generic.clear();
}
