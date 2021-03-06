
#include <utility>

#include <catch.hpp>

#include "../igraphpp/igraph.hpp"

#define N 5
#define M 7

TEST_CASE("Vector", "[Vector]") {
  using igraph::Vector;

  Vector vempty;
  CHECK(vempty.size() == 0);
  CHECK(vempty.empty());

  Vector vinit(10);
  CHECK(vinit.size() == 10);

  double a[N] = {1, 2, 3, 4, 5};
  double b[M] = {5, 5, 1, 4, 2, 3, 7};

  Vector va(a, N);
  Vector vb(a, N);

  Vector lst({1, 10});
  CHECK(lst.size() == 2);

  Vector seq(1, 5);
  CHECK(seq.size() == 5);
  CHECK(va == vb);

  CHECK(va.size() == 5);
  va.null();
  CHECK(va[2] == 0);
  va[1] = 4;
  CHECK(va.at(1) == 4);
  va.at(1) = 3;
  CHECK(va[1] == 3);
  CHECK(a[2] == 3);
  va.clear();
  CHECK(va.empty());

  va = Vector(a, N);
  va = vb + 2;
  va = va + 1;
  va = va - 1;
  va = vb * 2;
  vb = va / 2;
  vb = va;
  va = va;
  vb.append(va);
  CHECK(va.size() == N);
  CHECK(vb.size() == N + N);
  va.update(vb);
  vb.clear();
  CHECK(va.size() == N + N);
  CHECK(vb.empty());
  va.reverse();

  va = Vector(5);
  CHECK(va.size() == 5);
  va.reserve(10);
  va.resize(10);
  CHECK(va.size() == 10);
  CHECK(va.size() == N + N);
  va.push_back(-3);
  CHECK(va.size() == N + N + 1);
  CHECK(va.pop_back() == -3);
  CHECK(va.size() == N + N);
  va.push_back(1);
  va.remove(2);
  CHECK(va.size() == N + N);

  // vb.resize_min(); // This one never worked for me.

  va = Vector(a, N);
  vb = Vector(b, M);
  va.sort();
  vb.sort();

  Vector difference(va.difference_sorted(vb));
  CHECK(difference.size() == 0);
  difference = vb.difference_sorted(va);
  CHECK(difference.size() == 1);
  CHECK(difference[0] == 7);

  Vector intersect(va.intersect_sorted(vb));
  CHECK(intersect.size() == 5);
  CHECK(intersect[0] == 1);
  CHECK(intersect[1] == 2);
  CHECK(intersect[2] == 3);
  CHECK(intersect[3] == 4);
  CHECK(intersect[4] == 5);

  difference = std::move(intersect);
  Vector x(std::move(difference));
  Vector s(std::move(x));

  Vector rep = Vector::Repeat(5, 4);
  CHECK(rep == Vector({5, 5, 5, 5}));
}
