#include "../igraphpp/igraph.hpp"

#include <igraph.h>

namespace igraph {

/* Constructors and Destructors */
VectorView::~VectorView() {
  // This |vector_| must not be released from the memory.
}
VectorView::VectorView(const double *data, long int length) {
  igraph_vector_view(ptr(), data, length);
}

/* Initializing elements */
void VectorView::null() noexcept { igraph_vector_null(ptr()); }
void VectorView::fill(double value) noexcept {
  igraph_vector_fill(ptr(), value);
}

/* Accessing elements */
double &VectorView::operator[](long int i) const noexcept {
  return VECTOR(vector_)[i];
}
double &VectorView::at(long int i) const noexcept { return VECTOR(vector_)[i]; }
double &VectorView::head() const noexcept { return this->at(0); }
double &VectorView::tail() const noexcept { return this->at(size() - 1); }

/* Copying vectors */
void VectorView::copy_to(double *output) noexcept {
  igraph_vector_copy_to(ptr(), output);
}

/* Exchanging elements */
void VectorView::swap(long int i, long int j) noexcept {
  SafeCall(igraph_vector_swap_elements(ptr(), i, j));
}
void VectorView::swap(VectorView &b) {
  SafeCall(igraph_vector_swap(ptr(), b.ptr()));
}
void VectorView::shuffle() noexcept { SafeCall(igraph_vector_shuffle(ptr())); }
void VectorView::reverse() noexcept { SafeCall(igraph_vector_reverse(ptr())); }

/* Vector operations */
VectorView &VectorView::operator+=(double scalar) noexcept {
  igraph_vector_add_constant(ptr(), scalar);
  return *this;
}
VectorView &VectorView::operator-=(double scalar) noexcept {
  *this += -scalar;
  return *this;
}
VectorView &VectorView::operator*=(double scalar) noexcept {
  igraph_vector_scale(ptr(), scalar);
  return *this;
}
VectorView &VectorView::operator/=(double scalar) noexcept {
  return this->operator*=(1.0 / scalar);
}
VectorView &VectorView::operator+=(const VectorView &b) {
  SafeCall(igraph_vector_add(ptr(), b.ptr()));
  return *this;
}
VectorView &VectorView::operator-=(const VectorView &b) {
  SafeCall(igraph_vector_sub(ptr(), b.ptr()));
  return *this;
}
VectorView &VectorView::operator*=(const VectorView &b) {
  SafeCall(igraph_vector_mul(ptr(), b.ptr()));
  return *this;
}
VectorView &VectorView::operator/=(const VectorView &b) {
  SafeCall(igraph_vector_div(ptr(), b.ptr()));
  return *this;
}

/* Vector comparisons */ bool VectorView::operator==(const VectorView &b) const
    noexcept {
  return igraph_vector_all_e(ptr(), b.ptr());
}
bool VectorView::operator!=(const VectorView &b) const noexcept {
  return !this->operator==(b);
}
bool VectorView::operator<(const VectorView &b) const noexcept {
  return igraph_vector_all_l(ptr(), b.ptr());
}
bool VectorView::operator>(const VectorView &b) const noexcept {
  return igraph_vector_all_g(ptr(), b.ptr());
}
bool VectorView::operator<=(const VectorView &b) const noexcept {
  return igraph_vector_all_le(ptr(), b.ptr());
}
bool VectorView::operator>=(const VectorView &b) const noexcept {
  return igraph_vector_all_ge(ptr(), b.ptr());
}

/* Finding minimum and maximum */
double VectorView::min() const noexcept { return igraph_vector_min(ptr()); }
double VectorView::max() const noexcept { return igraph_vector_max(ptr()); }
long int VectorView::which_min() const noexcept {
  return igraph_vector_which_min(ptr());
}
long int VectorView::which_max() const noexcept {
  return igraph_vector_which_max(ptr());
}
void VectorView::minmax(double *min, double *max) const {
  SafeCall(igraph_vector_minmax(ptr(), min, max));
}
void VectorView::which_minmax(long int *min, long int *max) const {
  SafeCall(igraph_vector_which_minmax(ptr(), min, max));
}
/* Vector properties */
bool VectorView::empty() const noexcept { return igraph_vector_empty(ptr()); }
long int VectorView::size() const noexcept { return igraph_vector_size(ptr()); }
long int VectorView::capacity() const noexcept {
  return igraph_vector_capacity(ptr());
}
double VectorView::sum() const noexcept { return igraph_vector_sum(ptr()); }
double VectorView::prod() const noexcept { return igraph_vector_prod(ptr()); }
bool VectorView::isininterval(double low, double high) const noexcept {
  return igraph_vector_isininterval(ptr(), low, high);
}
double VectorView::maxdifference(const VectorView &b) const noexcept {
  return igraph_vector_maxdifference(ptr(), b.ptr());
}
/* Searching for elements */
bool VectorView::contains(double value) const noexcept {
  return igraph_vector_contains(ptr(), value);
}
long int VectorView::search(double value, long int from) const noexcept {
  long int result;
  bool contains = igraph_vector_search(ptr(), from, value, &result);
  if (contains) return result;
  return -1;
}
long int VectorView::binsearch(double value) const noexcept {
  long int result;
  bool contains = igraph_vector_binsearch(ptr(), value, &result);
  if (contains) return result;
  return -1;
}

/* Sorting */
void VectorView::sort() noexcept { igraph_vector_sort(ptr()); }

Vector::~Vector() {
  if (VECTOR(*ptr()) != NULL) igraph_vector_destroy(ptr());
}
Vector::Vector(long int size) : VectorView() {
  SafeCall(igraph_vector_init(ptr(), size));
}
Vector::Vector(const double *data, long int length) : VectorView() {
  igraph_real_t *dd = const_cast<igraph_real_t *>(data);
  SafeCall(igraph_vector_init_copy(ptr(), dd, length));
}
Vector::Vector(double from, double to) : VectorView() {
  SafeCall(igraph_vector_init_seq(ptr(), from, to));
}
Vector::Vector(std::initializer_list<double> list) {
  SafeCall(igraph_vector_init(ptr(), 0));
  reserve(static_cast<long int>(list.size()));
  auto begin = list.begin();
  auto end = list.end();
  for (; begin != end; ++begin) {
    push_back(static_cast<double>(*begin));
  }
}
Vector::Vector(const Vector &other)
    : Vector(static_cast<const VectorView &>(other)) {}
Vector::Vector(const VectorView &other) : VectorView(other.is_none()) {
  if (!other.is_none()) SafeCall(igraph_vector_copy(ptr(), other.ptr()));
}
Vector::Vector(Vector &&other) : VectorView(other.is_none()) {
  *ptr() = *other.ptr();
  other.disown();
}
Vector &Vector::operator=(Vector &&other) & {
  if (owner()) igraph_vector_destroy(ptr());
  *ptr() = *other.ptr();
  other.disown();
  return *this;
}

/* Vector operations */
Vector Vector::operator+(double scalar) const noexcept {
  Vector a(*this);
  igraph_vector_add_constant(a.ptr(), scalar);
  return a;
}
Vector Vector::operator-(double scalar) const noexcept {
  return operator+(-scalar);
}
Vector Vector::operator*(double scalar) const noexcept {
  Vector a(*this);
  igraph_vector_scale(a.ptr(), scalar);
  return a;
}
Vector Vector::operator/(double scalar) const noexcept {
  return this->operator*(1.0 / scalar);
}
Vector Vector::operator+(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_add(a.ptr(), b.ptr()));
  return a;
}
Vector Vector::operator-(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_sub(a.ptr(), b.ptr()));
  return a;
}
Vector Vector::operator*(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_mul(a.ptr(), b.ptr()));
  return a;
}
Vector Vector::operator/(const VectorView &b) const {
  Vector a(*this);
  SafeCall(igraph_vector_div(a.ptr(), b.ptr()));
  return a;
}

/* Copying vectors */
Vector &Vector::operator=(const VectorView &other) & {
  if (this == &other) return *this;
  if (owner() && !is_none()) this->~Vector();
  SafeCall(igraph_vector_copy(ptr(), other.ptr()));
  return *this;
}
Vector &Vector::operator=(const Vector &other) & {
  if (this == &other) return *this;
  if (owner() && !is_none()) this->~Vector();
  SafeCall(igraph_vector_copy(ptr(), other.ptr()));
  return *this;
}
void Vector::append(const VectorView &other) {
  SafeCall(igraph_vector_append(ptr(), other.ptr()));
}
void Vector::update(const VectorView &other) noexcept {
  igraph_vector_update(ptr(), other.ptr());
}

/* Resizing operations */
void Vector::clear() noexcept { igraph_vector_clear(ptr()); }
void Vector::reserve(long int newsize) {
  SafeCall(igraph_vector_reserve(ptr(), newsize));
}
void Vector::resize(long int newsize) {
  SafeCall(igraph_vector_resize(ptr(), newsize));
}
void Vector::resize_min() { SafeCall(igraph_vector_resize_min(ptr())); }
void Vector::push_back(double value) {
  SafeCall(igraph_vector_push_back(ptr(), value));
}
double Vector::pop_back() noexcept { return igraph_vector_pop_back(ptr()); }
void Vector::insert(long int pos, double value) {
  SafeCall(igraph_vector_insert(ptr(), pos, value));
}
void Vector::remove(long int pos) noexcept { igraph_vector_remove(ptr(), pos); }
void Vector::remove(long int from, long int to) noexcept {
  igraph_vector_remove_section(ptr(), from, to);
}

/* Set operations on sorted vectors */
Vector Vector::intersect_sorted(const VectorView &other) const {
  Vector result;
  SafeCall(igraph_vector_intersect_sorted(ptr(), other.ptr(), result.ptr()));
  return result;
}
Vector Vector::difference_sorted(const VectorView &other) const {
  Vector result;
  SafeCall(igraph_vector_difference_sorted(ptr(), other.ptr(), result.ptr()));
  return result;
}

Vector Vector::Repeat(double value, long int times) {
  Vector vector(times);
  vector.fill(value);
  return vector;
}

Vector::Vector(const igraph_vector_t &vector) : VectorView() {
  SafeCall(igraph_vector_copy(ptr(), &vector));
}

}  // namespace igraph
