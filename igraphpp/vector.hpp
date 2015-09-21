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
  VectorView(const double *data, long int length);

  ~VectorView();

  double &operator[](long int i) const noexcept;

  bool operator==(const VectorView &b) const noexcept;

  bool operator<(const VectorView &b) const noexcept;

  bool operator>(const VectorView &b) const noexcept;

  bool operator<=(const VectorView &b) const noexcept;

  bool operator>=(const VectorView &b) const noexcept;

  double min() const noexcept;

  double max() const noexcept;

  double sum() const noexcept;

  double prod() const noexcept;

  long int which_min() const noexcept;

  long int which_max() const noexcept;

  void minmax(double *min, double *max) const;

  void which_minmax(long int *min, long int *max) const;

  bool IsInInterval(double low, double high) const noexcept;

  double MaxDifference(const VectorView &b) const noexcept;

  bool Contains(double value) const noexcept;

  long int Search(double value, long int from = 0) const noexcept;

  long int BinarySearch(double value) const noexcept;

  void Sort() noexcept;

  void SetNull() noexcept;

  void Fill(double value) noexcept;

  bool empty() const noexcept;

  long int size() const noexcept;

  long int capacity() const noexcept;

  const igraph_vector_t *ptr() const;

  void CopyTo(double *output) noexcept;

protected:
  VectorView() = default;

  igraph_vector_t vector_;
};

class Vector : public VectorView {
public:
  Vector(int long size = 0);

  Vector(const double *data, long int length);

  Vector(double from, double to);

  ~Vector();

  Vector(const Vector &other);

  Vector(Vector &&other);

  Vector &operator=(const Vector &other);

  Vector &operator=(Vector &&other);

  Vector operator+(double scalar) noexcept;

  Vector operator-(double scalar) noexcept;

  Vector operator*(double scalar) noexcept;

  Vector &operator+=(double scalar) noexcept;

  Vector &operator-=(double scalar) noexcept;

  Vector &operator*=(double scalar) noexcept;

  Vector operator*(const VectorView &b);

  Vector operator/(const VectorView &b);

  Vector &operator+=(const VectorView &b);

  Vector &operator-=(const VectorView &b);

  Vector &operator*=(const VectorView &b);

  Vector &operator/=(const VectorView &b);

  igraph_vector_t *ptr();

  const igraph_vector_t *ptr() const;

  void Clear() noexcept;

  void Reserve(long int size);

  void Resize(long int size);

  void ResizeMin();

  void PushBack(double value);

  double PopBack() noexcept;

  void Insert(long int pos, double value);

  void Remove(long int pos) noexcept;

  void Remove(long int from, long int to) noexcept;

  void CopyTo(double *output) noexcept;

  void Update(const Vector &input) noexcept;

  void Append(const Vector &input);

  void Swap(long int i, long int j) noexcept;

  void Swap(Vector &b);

  void Reverse() noexcept;

  void Shuffle() noexcept;

protected:
  explicit Vector(const igraph_vector_t &vector);
};

} // namespace igraph

#include "./vector_impl.hpp"

#endif // IGRAPHPP_VECTOR_HPP_
