#ifndef IGRAPHPP_VECTOR_HPP_
#define IGRAPHPP_VECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <iostream> // TODO remove

#include <igraph.h>

#include "./exception.hpp"

namespace igraph {

class Vector;

class VectorView {
public:
  /* Constructors and Destructors */
  ~VectorView() {
    // This |vector_| must not be released from the memory.
  }
  VectorView(const double *data, long int length) {
    igraph_vector_view(ptr(), data, length);
  }
  VectorView &operator=(VectorView &other) = delete;
  VectorView &operator=(VectorView &&other) = delete;

  /* Initializing elements */
  void null() noexcept { igraph_vector_null(ptr()); }
  void fill(double value) noexcept { igraph_vector_fill(ptr(), value); }

  /* Accessing elements */
  double &operator[](long int i) const noexcept { return VECTOR(vector_)[i]; }
  double &at(long int i) const noexcept { return VECTOR(vector_)[i]; }

  /* Copying vectors */
  void copy_to(double *output) noexcept {
    igraph_vector_copy_to(ptr(), output);
  }

  /* Exchanging elements */
  void swap(long int i, long int j) noexcept {
    SafeCall(igraph_vector_swap_elements(ptr(), i, j));
  }
  void swap(VectorView &b) { SafeCall(igraph_vector_swap(ptr(), b.ptr())); }
  void shuffle() noexcept { SafeCall(igraph_vector_shuffle(ptr())); }
  void reverse() noexcept { SafeCall(igraph_vector_reverse(ptr())); }

  /* Vector operations */
  VectorView &operator+=(double scalar) noexcept {
    igraph_vector_add_constant(ptr(), scalar);
    return *this;
  }
  VectorView &operator-=(double scalar) noexcept {
    *this += -scalar;
    return *this;
  }
  VectorView &operator*=(double scalar) noexcept {
    igraph_vector_scale(ptr(), scalar);
    return *this;
  }
  VectorView &operator/=(double scalar) noexcept {
    return this->operator*=(1.0 / scalar);
  }
  VectorView &operator+=(const VectorView &b) {
    SafeCall(igraph_vector_add(ptr(), b.ptr()));
    return *this;
  }
  VectorView &operator-=(const VectorView &b) {
    SafeCall(igraph_vector_sub(ptr(), b.ptr()));
    return *this;
  }
  VectorView &operator*=(const VectorView &b) {
    SafeCall(igraph_vector_mul(ptr(), b.ptr()));
    return *this;
  }
  VectorView &operator/=(const VectorView &b) {
    SafeCall(igraph_vector_div(ptr(), b.ptr()));
    return *this;
  }

  /* Vector comparisons */
  bool operator==(const VectorView &b) const noexcept {
    return igraph_vector_all_e(ptr(), b.ptr());
  }
  bool operator!=(const VectorView &b) const noexcept {
    return !this->operator==(b);
  }
  bool operator<(const VectorView &b) const noexcept {
    return igraph_vector_all_l(ptr(), b.ptr());
  }
  bool operator>(const VectorView &b) const noexcept {
    return igraph_vector_all_g(ptr(), b.ptr());
  }
  bool operator<=(const VectorView &b) const noexcept {
    return igraph_vector_all_le(ptr(), b.ptr());
  }
  bool operator>=(const VectorView &b) const noexcept {
    return igraph_vector_all_ge(ptr(), b.ptr());
  }

  /* Finding minimum and maximum */
  double min() const noexcept { return igraph_vector_min(ptr()); }
  double max() const noexcept { return igraph_vector_max(ptr()); }
  long int which_min() const noexcept { return igraph_vector_which_min(ptr()); }
  long int which_max() const noexcept { return igraph_vector_which_max(ptr()); }
  void minmax(double *min, double *max) const {
    SafeCall(igraph_vector_minmax(ptr(), min, max));
  }
  void which_minmax(long int *min, long int *max) const {
    SafeCall(igraph_vector_which_minmax(ptr(), min, max));
  }
  /* Vector properties */
  bool empty() const noexcept { return igraph_vector_empty(ptr()); }
  long int size() const noexcept { return igraph_vector_size(ptr()); }
  long int capacity() const noexcept { return igraph_vector_capacity(ptr()); }
  double sum() const noexcept { return igraph_vector_sum(ptr()); }
  double prod() const noexcept { return igraph_vector_prod(ptr()); }
  bool isininterval(double low, double high) const noexcept {
    return igraph_vector_isininterval(ptr(), low, high);
  }
  double maxdifference(const VectorView &b) const noexcept {
    return igraph_vector_maxdifference(ptr(), b.ptr());
  }
  /* Searching for elements */
  bool contains(double value) const noexcept {
    return igraph_vector_contains(ptr(), value);
  }
  long int search(double value, long int from = 0) const noexcept {
    long int result;
    bool contains = igraph_vector_search(ptr(), from, value, &result);
    if (contains)
      return result;
    return -1;
  }
  long int binsearch(double value) const noexcept {
    long int result;
    bool contains = igraph_vector_binsearch(ptr(), value, &result);
    if (contains)
      return result;
    return -1;
  }

  /* Sorting */
  void sort() noexcept { igraph_vector_sort(ptr()); }

  /* Internal use */
  const igraph_vector_t *ptr() const { return &vector_; }
  igraph_vector_t *ptr() { return &vector_; }

protected:
  VectorView() = default;

private:
  igraph_vector_t vector_;

private:
  VectorView(const VectorView &) = default;
};

class Vector : public VectorView {
public:
  /* Constructors and Destructors */
  ~Vector() {
    if (VECTOR(*ptr()) != NULL)
      igraph_vector_destroy(ptr());
  }
  explicit Vector(int long size = 0) : VectorView() {
    SafeCall(igraph_vector_init(ptr(), size));
  }
  Vector(const double *data, long int length) : VectorView() {
    igraph_real_t *dd = const_cast<igraph_real_t *>(data);
    SafeCall(igraph_vector_init_copy(ptr(), dd, length));
  }
  Vector(double from, double to) : VectorView() {
    SafeCall(igraph_vector_init_seq(ptr(), from, to));
  }
  Vector(const Vector &other) : VectorView() {
    SafeCall(igraph_vector_copy(ptr(), other.ptr()));
  }
  Vector(const VectorView &other) : VectorView() {
    SafeCall(igraph_vector_copy(ptr(), other.ptr()));
  }
  Vector(Vector &&other) : VectorView() {
    *ptr() = *other.ptr();
    VECTOR(*other.ptr()) = NULL;
  }
  Vector &operator=(Vector &&other) {
    this->~Vector();
    *ptr() = *other.ptr();
    VECTOR(*other.ptr()) = NULL;
    return *this;
  }

  /* Vector operations */
  Vector operator+(double scalar) const noexcept {
    Vector a(*this);
    igraph_vector_add_constant(a.ptr(), scalar);
    return a;
  }
  Vector operator-(double scalar) const noexcept { return operator+(-scalar); }
  Vector operator*(double scalar) const noexcept {
    Vector a(*this);
    igraph_vector_scale(a.ptr(), scalar);
    return a;
  }
  Vector operator/(double scalar) const noexcept {
    return this->operator*(1.0 / scalar);
  }
  Vector operator+(const VectorView &b) const {
    Vector a(*this);
    SafeCall(igraph_vector_add(a.ptr(), b.ptr()));
    return a;
  }
  Vector operator-(const VectorView &b) const {
    Vector a(*this);
    SafeCall(igraph_vector_sub(a.ptr(), b.ptr()));
    return a;
  }
  Vector operator*(const VectorView &b) const {
    Vector a(*this);
    SafeCall(igraph_vector_mul(a.ptr(), b.ptr()));
    return a;
  }
  Vector operator/(const VectorView &b) const {
    Vector a(*this);
    SafeCall(igraph_vector_div(a.ptr(), b.ptr()));
    return a;
  }

  /* Copying vectors */
  Vector &operator=(const VectorView &other) {
    igraph_vector_update(ptr(), other.ptr());
    return *this;
  }
  Vector &operator=(const Vector &other) {
    if (this != &other)
      igraph_vector_update(ptr(), other.ptr());
    return *this;
  }
  void append(const VectorView &other) {
    SafeCall(igraph_vector_append(ptr(), other.ptr()));
  }
  void update(const VectorView &other) noexcept {
    igraph_vector_update(ptr(), other.ptr());
  }

  /* Resizing operations */
  void clear() noexcept { igraph_vector_clear(ptr()); }
  void reserve(long int newsize) {
    SafeCall(igraph_vector_reserve(ptr(), newsize));
  }
  void resize(long int newsize) {
    SafeCall(igraph_vector_resize(ptr(), newsize));
  }
  void resize_min() { SafeCall(igraph_vector_resize_min(ptr())); }
  void push_back(double value) {
    SafeCall(igraph_vector_push_back(ptr(), value));
  }
  double pop_back() noexcept { return igraph_vector_pop_back(ptr()); }
  void insert(long int pos, double value) {
    SafeCall(igraph_vector_insert(ptr(), pos, value));
  }
  void remove(long int pos) noexcept { igraph_vector_remove(ptr(), pos); }
  void remove(long int from, long int to) noexcept {
    igraph_vector_remove_section(ptr(), from, to);
  }

  /* Set operations on sorted vectors */
  Vector intersect_sorted(const VectorView &other) const {
    Vector result;
    SafeCall(igraph_vector_intersect_sorted(ptr(), other.ptr(), result.ptr()));
    return result;
  }
  Vector difference_sorted(const VectorView &other) const {
    Vector result;
    SafeCall(igraph_vector_difference_sorted(ptr(), other.ptr(), result.ptr()));
    return result;
  }

private:
  explicit Vector(const igraph_vector_t &vector) : VectorView() {
    SafeCall(igraph_vector_copy(ptr(), &vector));
  }
};

} // namespace igraph

#endif // IGRAPHPP_VECTOR_HPP_
