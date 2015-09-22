#ifndef IGRAPHPP_VECTOR_HPP_
#define IGRAPHPP_VECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class VectorView {
public:
  /* Constructors and Destructors */
  ~VectorView();
  VectorView(const double *data, long int length);

  /* Initializing elements */
  void null() noexcept;
  void fill(double value) noexcept;

  /* Accessing elements */
  double &operator[](long int i) const noexcept;
  double &at(long int i) const noexcept;

  /* Copying vectors */
  void copy_to(double *output) noexcept;
  void update(const VectorView &input) noexcept;

  /* Exchanging elements */
  void swap(long int i, long int j) noexcept;
  void swap(VectorView &b);
  void shuffle() noexcept;
  void reverse() noexcept;

  /* Vector operations */
  VectorView operator+(double scalar) noexcept;
  VectorView operator-(double scalar) noexcept;
  VectorView operator*(double scalar) noexcept;
  VectorView &operator+=(double scalar) noexcept;
  VectorView &operator-=(double scalar) noexcept;
  VectorView &operator*=(double scalar) noexcept;
  VectorView operator*(const VectorView &b);
  VectorView operator/(const VectorView &b);
  VectorView &operator+=(const VectorView &b);
  VectorView &operator-=(const VectorView &b);
  VectorView &operator*=(const VectorView &b);
  VectorView &operator/=(const VectorView &b);

  /* Vector comparisons */
  bool operator==(const VectorView &b) const noexcept;
  bool operator<(const VectorView &b) const noexcept;
  bool operator>(const VectorView &b) const noexcept;
  bool operator<=(const VectorView &b) const noexcept;
  bool operator>=(const VectorView &b) const noexcept;

  /* Finding minimum and maximum */
  double min() const noexcept;
  double max() const noexcept;
  long int which_min() const noexcept;
  long int which_max() const noexcept;
  void minmax(double *min, double *max) const;
  void which_minmax(long int *min, long int *max) const;

  /* Vector properties */
  bool empty() const noexcept;
  long int size() const noexcept;
  long int capacity() const noexcept;
  double sum() const noexcept;
  double prod() const noexcept;
  bool isininterval(double low, double high) const noexcept;
  double maxdifference(const VectorView &b) const noexcept;

  /* Searching for elements */
  bool contains(double value) const noexcept;
  long int search(double value, long int from = 0) const noexcept;
  long int binsearch(double value) const noexcept;

  /* Sorting */
  void sort() noexcept;

  /* Set operations on sorted vectors */
  VectorView intersect_sorted(const VectorView &v2);
  VectorView difference_sorted(const VectorView &v2);

  /* Internal use */
  const igraph_vector_t *ptr() const;
  igraph_vector_t *ptr();

protected:
  VectorView() = default;

  igraph_vector_t vector_;
};

class Vector : public VectorView {
public:
  /* Constructors and Destructors */
  ~Vector();
  Vector(int long size = 0);
  Vector(const double *data, long int length);
  Vector(double from, double to);
  Vector(const Vector &other);
  Vector(Vector &&other);

  Vector &operator=(Vector &&other);

  /* Copying vectors */
  Vector &operator=(const Vector &other);
  void append(const Vector &input);

  /* Resizing operations */
  void clear() noexcept;
  void reserve(long int size);
  void resize(long int size);
  void resize_min();
  void push_back(double value);
  double pop_back() noexcept;
  void insert(long int pos, double value);
  void remove(long int pos) noexcept;
  void remove(long int from, long int to) noexcept;

protected:
  explicit Vector(const igraph_vector_t &vector);
};

} // namespace igraph

#include "./vector_impl.hpp"

#endif // IGRAPHPP_VECTOR_HPP_
