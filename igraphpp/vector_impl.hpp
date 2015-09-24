#ifndef IGRAPHPP_VECTOR_IMPL_HPP_
#define IGRAPHPP_VECTOR_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

/* Constructors and Destructors */
inline VectorView::~VectorView() {
  // This |vector_| must not be released from the memory.
}
inline VectorView::VectorView(const double *data, long int length) {
  igraph_vector_view(ptr(), data, length);
}

/* Initializing elements */
inline void VectorView::null() noexcept { igraph_vector_null(ptr()); }
inline void VectorView::fill(double value) noexcept {
  igraph_vector_fill(ptr(), value);
}

/* Accessing elements */
inline double &VectorView::operator[](long int i) const noexcept {
  return VECTOR(vector_)[i];
}
inline double &VectorView::at(long int i) const noexcept {
  return VECTOR(vector_)[i];
}
inline double &VectorView::head() const noexcept { return this->at(0); }
inline double &VectorView::tail() const noexcept {
  return this->at(size() - 1);
}

/* Copying vectors */
inline void VectorView::copy_to(double *output) noexcept {
  igraph_vector_copy_to(ptr(), output);
}

/* Exchanging elements */
inline void VectorView::swap(long int i, long int j) noexcept {
  SafeCall(igraph_vector_swap_elements(ptr(), i, j));
}
inline void VectorView::swap(VectorView &b) {
  SafeCall(igraph_vector_swap(ptr(), b.ptr()));
}
inline void VectorView::shuffle() noexcept {
  SafeCall(igraph_vector_shuffle(ptr()));
}
inline void VectorView::reverse() noexcept {
  SafeCall(igraph_vector_reverse(ptr()));
}

/* Vector operations */
inline VectorView &VectorView::operator+=(double scalar) noexcept {
  igraph_vector_add_constant(ptr(), scalar);
  return *this;
}
inline VectorView &VectorView::operator-=(double scalar) noexcept {
  *this += -scalar;
  return *this;
}
inline VectorView &VectorView::operator*=(double scalar) noexcept {
  igraph_vector_scale(ptr(), scalar);
  return *this;
}
inline VectorView &VectorView::operator/=(double scalar) noexcept {
  return this->operator*=(1.0 / scalar);
}
inline VectorView &VectorView::operator+=(const VectorView &b) {
  SafeCall(igraph_vector_add(ptr(), b.ptr()));
  return *this;
}
inline VectorView &VectorView::operator-=(const VectorView &b) {
  SafeCall(igraph_vector_sub(ptr(), b.ptr()));
  return *this;
}
inline VectorView &VectorView::operator*=(const VectorView &b) {
  SafeCall(igraph_vector_mul(ptr(), b.ptr()));
  return *this;
}
inline VectorView &VectorView::operator/=(const VectorView &b) {
  SafeCall(igraph_vector_div(ptr(), b.ptr()));
  return *this;
}

/* Vector comparisons */ inline bool VectorView::
operator==(const VectorView &b) const noexcept {
  return igraph_vector_all_e(ptr(), b.ptr());
}
inline bool VectorView::operator!=(const VectorView &b) const noexcept {
  return !this->operator==(b);
}
inline bool VectorView::operator<(const VectorView &b) const noexcept {
  return igraph_vector_all_l(ptr(), b.ptr());
}
inline bool VectorView::operator>(const VectorView &b) const noexcept {
  return igraph_vector_all_g(ptr(), b.ptr());
}
inline bool VectorView::operator<=(const VectorView &b) const noexcept {
  return igraph_vector_all_le(ptr(), b.ptr());
}
inline bool VectorView::operator>=(const VectorView &b) const noexcept {
  return igraph_vector_all_ge(ptr(), b.ptr());
}

/* Finding minimum and maximum */
inline double VectorView::min() const noexcept {
  return igraph_vector_min(ptr());
}
inline double VectorView::max() const noexcept {
  return igraph_vector_max(ptr());
}
inline long int VectorView::which_min() const noexcept {
  return igraph_vector_which_min(ptr());
}
inline long int VectorView::which_max() const noexcept {
  return igraph_vector_which_max(ptr());
}
inline void VectorView::minmax(double *min, double *max) const {
  SafeCall(igraph_vector_minmax(ptr(), min, max));
}
inline void VectorView::which_minmax(long int *min, long int *max) const {
  SafeCall(igraph_vector_which_minmax(ptr(), min, max));
}
/* Vector properties */
inline bool VectorView::empty() const noexcept {
  return igraph_vector_empty(ptr());
}
inline long int VectorView::size() const noexcept {
  return igraph_vector_size(ptr());
}
inline long int VectorView::capacity() const noexcept {
  return igraph_vector_capacity(ptr());
}
inline double VectorView::sum() const noexcept {
  return igraph_vector_sum(ptr());
}
inline double VectorView::prod() const noexcept {
  return igraph_vector_prod(ptr());
}
inline bool VectorView::isininterval(double low, double high) const noexcept {
  return igraph_vector_isininterval(ptr(), low, high);
}
inline double VectorView::maxdifference(const VectorView &b) const noexcept {
  return igraph_vector_maxdifference(ptr(), b.ptr());
}
/* Searching for elements */
inline bool VectorView::contains(double value) const noexcept {
  return igraph_vector_contains(ptr(), value);
}
inline long int VectorView::search(double value, long int from) const noexcept {
  long int result;
  bool contains = igraph_vector_search(ptr(), from, value, &result);
  if (contains)
    return result;
  return -1;
}
inline long int VectorView::binsearch(double value) const noexcept {
  long int result;
  bool contains = igraph_vector_binsearch(ptr(), value, &result);
  if (contains)
    return result;
  return -1;
}

/* Sorting */
inline void VectorView::sort() noexcept { igraph_vector_sort(ptr()); }

inline Vector::~Vector() {
  if (VECTOR(*ptr()) != NULL)
    igraph_vector_destroy(ptr());
}
inline Vector::Vector(int long size) : VectorView() {
  SafeCall(igraph_vector_init(ptr(), size));
}
inline Vector::Vector(const double *data, long int length) : VectorView() {
  igraph_real_t *dd = const_cast<igraph_real_t *>(data);
  SafeCall(igraph_vector_init_copy(ptr(), dd, length));
}
inline Vector::Vector(double from, double to) : VectorView() {
  SafeCall(igraph_vector_init_seq(ptr(), from, to));
}
inline Vector::Vector(const Vector &other)
    : Vector(static_cast<const VectorView &>(other)) {}
inline Vector::Vector(const VectorView &other) : VectorView() {
  SafeCall(igraph_vector_copy(ptr(), other.ptr()));
}
inline Vector::Vector(Vector &&other) : VectorView() {
  *ptr() = *other.ptr();
  VECTOR(*other.ptr()) = NULL;
}
inline Vector &Vector::operator=(Vector &&other) {
  this->~Vector();
  *ptr() = *other.ptr();
  VECTOR(*other.ptr()) = NULL;
  return *this;
}

/* Vector operations */
inline Vector Vector::operator+(double scalar) const noexcept {
  Vector a(*this);
  igraph_vector_add_constant(a.ptr(), scalar);
  return a;
}
inline Vector Vector::operator-(double scalar) const noexcept {
  return operator+(-scalar);
}
inline Vector Vector::operator*(double scalar) const noexcept {
  Vector a(*this);
  igraph_vector_scale(a.ptr(), scalar);
  return a;
}
inline Vector Vector::operator/(double scalar) const noexcept {
  return this->operator*(1.0 / scalar);
}
inline Vector Vector::operator+(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_add(a.ptr(), b.ptr()));
  return a;
}
inline Vector Vector::operator-(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_sub(a.ptr(), b.ptr()));
  return a;
}
inline Vector Vector::operator*(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_mul(a.ptr(), b.ptr()));
  return a;
}
inline Vector Vector::operator/(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_div(a.ptr(), b.ptr()));
  return a;
}

/* Copying vectors */
inline Vector &Vector::operator=(const VectorView &other) {
  igraph_vector_update(ptr(), other.ptr());
  return *this;
}
inline Vector &Vector::operator=(const Vector &other) {
  if (this != &other)
    igraph_vector_update(ptr(), other.ptr());
  return *this;
}
inline void Vector::append(const VectorView &other) {
  SafeCall(igraph_vector_append(ptr(), other.ptr()));
}
inline void Vector::update(const VectorView &other) noexcept {
  igraph_vector_update(ptr(), other.ptr());
}

/* Resizing operations */
inline void Vector::clear() noexcept { igraph_vector_clear(ptr()); }
inline void Vector::reserve(long int newsize) {
  SafeCall(igraph_vector_reserve(ptr(), newsize));
}
inline void Vector::resize(long int newsize) {
  SafeCall(igraph_vector_resize(ptr(), newsize));
}
inline void Vector::resize_min() { SafeCall(igraph_vector_resize_min(ptr())); }
inline void Vector::push_back(double value) {
  SafeCall(igraph_vector_push_back(ptr(), value));
}
inline double Vector::pop_back() noexcept {
  return igraph_vector_pop_back(ptr());
}
inline void Vector::insert(long int pos, double value) {
  SafeCall(igraph_vector_insert(ptr(), pos, value));
}
inline void Vector::remove(long int pos) noexcept {
  igraph_vector_remove(ptr(), pos);
}
inline void Vector::remove(long int from, long int to) noexcept {
  igraph_vector_remove_section(ptr(), from, to);
}

/* Set operations on sorted vectors */
inline Vector Vector::intersect_sorted(const VectorView &other) const {
  Vector result;
  SafeCall(igraph_vector_intersect_sorted(ptr(), other.ptr(), result.ptr()));
  return result;
}
inline Vector Vector::difference_sorted(const VectorView &other) const {
  Vector result;
  SafeCall(igraph_vector_difference_sorted(ptr(), other.ptr(), result.ptr()));
  return result;
}

inline Vector::Vector(const igraph_vector_t &vector) : VectorView() {
  SafeCall(igraph_vector_copy(ptr(), &vector));
}

} // namespace igraph

#endif // IGRAPHPP_VECTOR_IMPL_HPP_
