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
  VectorView(const double* data, long int length);

  ~VectorView();

  long int size() const noexcept;

  const igraph_vector_t* ptr() const;

private:
  igraph_vector_t vector_;
};

class Vector {
public:
  Vector(int long size);

  Vector(const double *data, long int length);

  Vector(double from, double to);

  ~Vector();

  Vector(const Vector &other);

  Vector(Vector &&other);

  Vector &operator=(const Vector &other);

  Vector &operator=(const Vector &&other);

  double operator[](long int i) const noexcept;

  void SetNull() noexcept;

  void Fill(double value) noexcept;

  void CopyTo(double *output) noexcept;

  void Update(const Vector& input) noexcept;

  void Append(const Vector& input);

  long int size() const noexcept;

protected:
  Vector(const igraph_vector_t &vector);

private:
  igraph_vector_t vector_;
};

} // namespace igraph

#include "./vector_impl.hpp"

#endif // IGRAPHPP_VECTOR_HPP_
