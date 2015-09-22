
#include <catch.hpp>
#include "../igraphpp/igraph.hpp"

#define N 42

TEST_CASE("VectorView", "[VectorView]") {
  using igraph::VectorView;

  double a[N];
  double b[N];

  VectorView va(a, N);
  VectorView vb(b, N);

  va.fill(1.0);
  vb.null();

  CHECK(va > vb);
  CHECK(vb < va);

  vb += 1.0;
  CHECK(va == vb);
  CHECK(va >= vb);

  va += 2.0;
  CHECK(va > vb);
  vb.update(va);

  va = vb * 4.0;
  va[0] = 10.0;
  CHECK(va[0] == Approx(10.0));
  CHECK(va.at(0) == Approx(10.0));
  CHECK(vb[0] != Approx(10.0));
  
  va[4] = -4.0;
  va.swap(0, 4);
  CHECK(va.at(0) == Approx(-4.0));
  CHECK(va[4] == Approx(10.0));

  va.swap(vb);
  CHECK(vb.at(0) == Approx(-4.0));
  CHECK(vb[4] == Approx(10.0));
  CHECK_FALSE(va.at(0) == Approx(-4.0));
  CHECK_FALSE(va[4] == Approx(10.0));

}
