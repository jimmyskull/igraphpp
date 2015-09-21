#ifndef IGRAPH_VECTOR_IMPL_HPP_
#define IGRAPH_VECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

VectorView::VectorView(const double* data, long int length) {
  igraph_vector_view(&vector_, data, length);
}

VectorView::~VectorView() {
  // This |vector_| must not be released from the memory.
}

long int VectorView::size() const noexcept {
  return igraph_vector_size(&vector_);
}

const igraph_vector_t* VectorView::ptr() const {
  return &vector_;
}

Vector::Vector(const igraph_vector_t &vector) {
  SafeCall(igraph_vector_copy(&vector_, &vector));
}

Vector::Vector(int long size) { SafeCall(igraph_vector_init(&vector_, size)); }

Vector::Vector(const double *data, long int length) {
  igraph_real_t *dd = const_cast<igraph_real_t *>(data);
  SafeCall(igraph_vector_init_copy(&vector_, dd, length));
}

Vector::Vector(double from, double to) {
  SafeCall(igraph_vector_init_seq(&vector_, from, to));
}

Vector::~Vector() {
  // TODO: Note that vectors created by igraph_vector_view() are special, you
  // mustn't call igraph_vector_destroy() on these.
  igraph_vector_destroy(&vector_);
}

Vector::Vector(const Vector &other) {
  SafeCall(igraph_vector_copy(&vector_, &other.vector_));
}

Vector::Vector(Vector &&other) { vector_ = std::move(other.vector_); }

Vector &Vector::operator=(const Vector &other) {
  if (this == &other)
    return *this;
  igraph_vector_destroy(&vector_);
  SafeCall(igraph_vector_copy(&vector_, &other.vector_));
  return *this;
}

Vector &Vector::operator=(const Vector &&other) {
  vector_ = std::move(other.vector_);
  return *this;
}

double Vector::operator[](long int i) const noexcept {
  return VECTOR(vector_)[i];
}

void Vector::SetNull() noexcept {
  igraph_vector_null(&vector_);
}

void Vector::Fill(double value) noexcept {
  igraph_vector_fill(&vector_, value);
}

void Vector::CopyTo(double *output) noexcept {
  igraph_vector_copy_to(&vector_, output);
}

void Vector::Update(const Vector& input) noexcept {
  igraph_vector_update(&vector_, &input.vector_);
}

void Vector::Append(const Vector& input) {
  SafeCall(igraph_vector_append(&vector_, &input.vector_));
}

long int Vector::size() const noexcept {
  return igraph_vector_size(&vector_);
}


} // namespace igraph

#endif // IGRAPH_VECTOR_IMPL_HPP_
