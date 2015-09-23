
#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

#define N 5

TEST_CASE("VectorView", "[VectorView]") {
  using igraph::VectorView;
  using igraph::Vector;

  double a[N];
  double b[N];

  VectorView va(a, N);
  VectorView vb(b, N);

  CHECK_FALSE(va.empty());
  CHECK(va.size() == N);

  // va = 1 1 1 1 1
  // vb = 0 0 0 0 0
  va.fill(1.0);
  vb.null();

  CHECK_FALSE(va == vb);
  CHECK(va != vb);
  CHECK(va > vb);
  CHECK_FALSE(va < vb);
  CHECK(va >= vb);
  CHECK_FALSE(va <= vb);

  // va = 1 1 1 1 1
  // vb = 1 1 1 1 1
  vb += 1.0;
  CHECK(va == vb);
  CHECK_FALSE(va != vb);
  CHECK_FALSE(va > vb);
  CHECK_FALSE(va < vb);
  CHECK(va >= vb);
  CHECK(va <= vb);

  // va = 3 3 3 3 3
  // vb = 1 1 1 1 1
  va += 2.0;
  va /= 1.0;
  CHECK_FALSE(va == vb);
  CHECK(va != vb);
  CHECK(va > vb);
  CHECK_FALSE(va < vb);
  CHECK(va >= vb);
  CHECK_FALSE(va <= vb);

  // va = 3 3 3 3 3
  // vb = 3 3 3 3 3
  va.copy_to(b);
  CHECK(va == vb);
  CHECK_FALSE(va != vb);
  CHECK_FALSE(va > vb);
  CHECK_FALSE(va < vb);
  CHECK(va >= vb);
  CHECK(va <= vb);

  // va = 9 9 9 9 9
  // vb = 2 2 2 2 2
  va *= 3.0;
  vb -= 1.0;
  CHECK(va[0] == 9.0);
  CHECK(vb[0] == 2.0);

  // va = 9 9 9 9 9
  // vb = 1 2 2 2 2
  b[0] = 1.0;
  CHECK(va[0] == 9.0);
  CHECK(vb[0] == 1.0);
  CHECK_FALSE(va == vb);
  CHECK(va != vb);
  CHECK(va > vb);
  CHECK_FALSE(va < vb);
  CHECK(va >= vb);
  CHECK_FALSE(va <= vb);

  // va = 10 9 9 9 9
  // vb =  1 2 2 2 2
  va[0] = 10.0;
  CHECK(va[0] == 10.0);
  CHECK(va.at(0) == 10.0);
  CHECK(vb[0] != 10.0);

  // va = -4 9 9 9 10
  // vb =  1 2 2 2  2
  va[4] = -4.0;
  va.swap(0, 4);
  CHECK(va.at(0) == -4.0);
  CHECK(va[4] == 10.0);

  // va =  1 2 2 2  2
  // vb = -4 9 9 9 10
  va.swap(vb);
  CHECK(vb.at(0) == -4.0);
  CHECK(vb[4] == 10.0);
  CHECK_FALSE(va.at(0) == -4.0);
  CHECK_FALSE(va[4] == 10.0);
  CHECK(va.at(0) == 1.0);
  CHECK(va[4] == 2.0);

  CHECK(vb.sum() == 33);
  CHECK(vb.prod() == (-4 * 9 * 9 * 9 * 10.0));

  // va =  1 2 3 2  2
  // vb = -4 9 9 6 10
  va[2] = 3;
  vb[3] = 6;
  CHECK(va.isininterval(1, 3));
  CHECK(vb.isininterval(-4, 10));
  CHECK_FALSE(va.isininterval(2, 3));
  CHECK_FALSE(vb.isininterval(-4, 9));
  CHECK(vb.maxdifference(va) == 8);

  CHECK_FALSE(va.contains(0));
  CHECK(va.contains(3));
  CHECK(vb.contains(-4));

  // va = 1 2 2 3  4
  // vb = 1 2 6 9 10
  va.sort();
  vb.sort();
  vb[vb.binsearch(-4)] = 1;
  va[va.search(2)] = 4;
  vb[vb.binsearch(9)] = 2;
  va.sort();
  vb.sort();
  CHECK(vb.search(51) == -1);
  CHECK(vb.binsearch(51) == -1);

  CHECK(va <= vb);

  vb.reverse();
  CHECK(vb[0] == 10);
  CHECK(vb[N - 1] == 1);

  CHECK(vb.min() == 1);
  CHECK(vb.which_min() == N - 1);
  CHECK(vb.max() == 10);
  CHECK(vb.which_max() == 0);

  double min, max;
  long int wmin, wmax;
  va.minmax(&min, &max);
  va.which_minmax(&wmin, &wmax);
  CHECK(min == 1);
  CHECK(wmin == 0);
  CHECK(max == 4);
  CHECK(wmax == N - 1);

  // Vector intersect = va.intersect_sorted(vb);
  // CHECK(intersect.size() == 2);
  // CHECK(intersect[0] == 1);
  // CHECK(intersect[1] == 2);

  // Vector difference = va.difference_sorted(vb);
  // CHECK(difference.size() == 2);
  // CHECK(difference[0] == 3);
  // CHECK(difference[1] == 4);
}
