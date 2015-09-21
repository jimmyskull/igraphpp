#ifndef IGRAPH_VECTOR_IMPL_HPP_
#define IGRAPH_VECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

VectorView::VectorView(const double *data, long int length) {
  igraph_vector_view(&vector_, data, length);
}

VectorView::~VectorView() {
  // This |vector_| must not be released from the memory.
}

double &VectorView::operator[](long int i) const noexcept {
  return VECTOR(vector_)[i];
}

bool VectorView::operator==(const VectorView &b) const noexcept {
  return igraph_vector_all_e(&vector_, b.ptr());
}

bool VectorView::operator<(const VectorView &b) const noexcept {
  return igraph_vector_all_l(&vector_, b.ptr());
}

bool VectorView::operator>(const VectorView &b) const noexcept {
  return igraph_vector_all_g(&vector_, b.ptr());
}

bool VectorView::operator<=(const VectorView &b) const noexcept {
  return igraph_vector_all_le(&vector_, b.ptr());
}

bool VectorView::operator>=(const VectorView &b) const noexcept {
  return igraph_vector_all_ge(&vector_, b.ptr());
}

double VectorView::min() const noexcept { return igraph_vector_min(&vector_); }

double VectorView::max() const noexcept { return igraph_vector_max(&vector_); }

double VectorView::sum() const noexcept { return igraph_vector_sum(&vector_); }

double VectorView::prod() const noexcept {
  return igraph_vector_prod(&vector_);
}

long int VectorView::which_min() const noexcept {
  return igraph_vector_which_min(&vector_);
}

long int VectorView::which_max() const noexcept {
  return igraph_vector_which_max(&vector_);
}

void VectorView::minmax(double *min, double *max) const {
  SafeCall(igraph_vector_minmax(&vector_, min, max));
}

void VectorView::which_minmax(long int *min, long int *max) const {
  SafeCall(igraph_vector_which_minmax(&vector_, min, max));
}

bool VectorView::IsInInterval(double low, double high) const noexcept {
  return igraph_vector_isininterval(&vector_, low, high);
}

double VectorView::MaxDifference(const VectorView &b) const noexcept {
  return igraph_vector_maxdifference(&vector_, b.ptr());
}

bool VectorView::Contains(double value) const noexcept {
  return igraph_vector_contains(&vector_, value);
}

long int VectorView::Search(double value, long int from) const noexcept {
  long int result;
  bool contains = igraph_vector_search(&vector_, from, value, &result);
  if (contains)
    return result;
  return -1;
}

long int VectorView::BinarySearch(double value) const noexcept {
  long int result;
  bool contains = igraph_vector_binsearch(&vector_, value, &result);
  if (contains)
    return result;
  return -1;
}

void VectorView::Sort() noexcept { igraph_vector_sort(&vector_); }

void VectorView::SetNull() noexcept { igraph_vector_null(&vector_); }

void VectorView::Fill(double value) noexcept {
  igraph_vector_fill(&vector_, value);
}

bool VectorView::empty() const noexcept {
  return igraph_vector_empty(&vector_);
}

long int VectorView::size() const noexcept {
  return igraph_vector_size(&vector_);
}

long int VectorView::capacity() const noexcept {
  return igraph_vector_capacity(&vector_);
}

const igraph_vector_t *VectorView::ptr() const { return &vector_; }

void VectorView::CopyTo(double *output) noexcept {
  igraph_vector_copy_to(&vector_, output);
}

Vector::Vector(const igraph_vector_t &vector) : VectorView() {
  SafeCall(igraph_vector_copy(&vector_, &vector));
}

Vector::Vector(int long size) : VectorView() {
  SafeCall(igraph_vector_init(&vector_, size));
}

Vector::Vector(const double *data, long int length) : VectorView() {
  igraph_real_t *dd = const_cast<igraph_real_t *>(data);
  SafeCall(igraph_vector_init_copy(&vector_, dd, length));
}

Vector::Vector(double from, double to) : VectorView() {
  SafeCall(igraph_vector_init_seq(&vector_, from, to));
}

Vector::~Vector() {
  if (VECTOR(vector_) != NULL)
    igraph_vector_destroy(&vector_);
}

Vector::Vector(const Vector &other) : VectorView() {
  SafeCall(igraph_vector_copy(&vector_, &other.vector_));
}

Vector::Vector(Vector &&other) : VectorView() {
  vector_ = other.vector_;
  VECTOR(other.vector_) = NULL;
}

Vector &Vector::operator=(const Vector &other) {
  if (this == &other)
    return *this;
  igraph_vector_destroy(&vector_);
  SafeCall(igraph_vector_copy(&vector_, &other.vector_));
  return *this;
}

Vector &Vector::operator=(Vector &&other) {
  vector_ = other.vector_;
  VECTOR(other.vector_) = NULL;
  return *this;
}

Vector Vector::operator+(double scalar) noexcept {
  Vector a(*this);
  igraph_vector_add_constant(a.ptr(), scalar);
  return a;
}

Vector Vector::operator-(double scalar) noexcept {
  Vector a(*this);
  igraph_vector_add_constant(a.ptr(), -scalar);
  return a;
}

Vector Vector::operator*(double scalar) noexcept {
  Vector a(*this);
  igraph_vector_scale(a.ptr(), scalar);
  return a;
}

Vector &Vector::operator+=(double scalar) noexcept {
  igraph_vector_add_constant(&vector_, scalar);
  return *this;
}

Vector &Vector::operator-=(double scalar) noexcept {
  igraph_vector_add_constant(&vector_, scalar);
  return *this;
}

Vector &Vector::operator*=(double scalar) noexcept {
  igraph_vector_scale(&vector_, scalar);
  return *this;
}

Vector Vector::operator*(const VectorView &b) {
  Vector a(*this);
  SafeCall(igraph_vector_mul(a.ptr(), b.ptr()));
  return a;
}

Vector Vector::operator/(const VectorView &b) {
  Vector a(*this);
  SafeCall(igraph_vector_div(&a.vector_, b.ptr()));
  return a;
}

Vector &Vector::operator+=(const VectorView &b) {
  SafeCall(igraph_vector_add(&vector_, b.ptr()));
  return *this;
}

Vector &Vector::operator-=(const VectorView &b) {
  SafeCall(igraph_vector_sub(&vector_, b.ptr()));
  return *this;
}

Vector &Vector::operator*=(const VectorView &b) {
  SafeCall(igraph_vector_mul(&vector_, b.ptr()));
  return *this;
}

Vector &Vector::operator/=(const VectorView &b) {
  SafeCall(igraph_vector_div(&vector_, b.ptr()));
  return *this;
}

igraph_vector_t *Vector::ptr() { return &vector_; }

const igraph_vector_t *Vector::ptr() const { return &vector_; }

void Vector::Clear() noexcept { igraph_vector_clear(&vector_); }

void Vector::Reserve(long int size) {
  SafeCall(igraph_vector_reserve(&vector_, size));
}

void Vector::Resize(long int size) {
  SafeCall(igraph_vector_resize(&vector_, size));
}

void Vector::ResizeMin() { SafeCall(igraph_vector_resize_min(&vector_)); }

void Vector::PushBack(double value) {
  SafeCall(igraph_vector_push_back(&vector_, value));
}

double Vector::PopBack() noexcept { return igraph_vector_pop_back(&vector_); }

void Vector::Insert(long int pos, double value) {
  SafeCall(igraph_vector_insert(&vector_, pos, value));
}

void Vector::Remove(long int pos) noexcept {
  igraph_vector_remove(&vector_, pos);
}

void Vector::Remove(long int from, long int to) noexcept {
  igraph_vector_remove_section(&vector_, from, to);
}

void Vector::CopyTo(double *output) noexcept {
  igraph_vector_copy_to(&vector_, output);
}

void Vector::Update(const Vector &input) noexcept {
  igraph_vector_update(&vector_, &input.vector_);
}

void Vector::Append(const Vector &input) {
  SafeCall(igraph_vector_append(&vector_, &input.vector_));
}

void Vector::Swap(long int i, long int j) noexcept {
  SafeCall(igraph_vector_swap_elements(&vector_, i, j));
}

void Vector::Swap(Vector &b) {
  SafeCall(igraph_vector_swap(&vector_, b.ptr()));
}

void Vector::Reverse() noexcept { SafeCall(igraph_vector_reverse(&vector_)); }

void Vector::Shuffle() noexcept { SafeCall(igraph_vector_shuffle(&vector_)); }

} // namespace igraph

#endif // IGRAPH_VECTOR_IMPL_HPP_
