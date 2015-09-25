#ifndef IGRAPHPP_MATRIX_IMPL_HPP_
#define IGRAPHPP_MATRIX_IMPL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include "./matrix.hpp"

#include <igraph.h>

namespace igraph {

/* Constructors and destructors */
inline Matrix::~Matrix() {
  if (owner())
    igraph_matrix_destroy(ptr());
}
inline Matrix::Matrix(long int nrow, long int ncol) {
  SafeCall(igraph_matrix_init(ptr(), nrow, ncol));
}
inline Matrix::Matrix(const Matrix &matrix) {
  SafeCall(igraph_matrix_copy(ptr(), matrix.ptr()));
}
inline Matrix::Matrix(Matrix &&matrix) {
  matrix_ = *matrix.ptr();
  matrix.disown();
}
inline Matrix &Matrix::operator=(const Matrix &matrix) {
  this->~Matrix();
  SafeCall(igraph_matrix_copy(ptr(), matrix.ptr()));
  return *this;
}
inline Matrix &Matrix::operator=(Matrix &&matrix) {
  this->~Matrix();
  *ptr() = *matrix.ptr();
  matrix.disown();
  return *this;
}

/* Initializing elements */
inline void Matrix::null() noexcept { igraph_matrix_null(ptr()); }
inline void Matrix::fill(double value) noexcept {
  igraph_matrix_fill(ptr(), value);
}

/* Copying matrices */
inline void Matrix::copy_to(double *to) noexcept {
  igraph_matrix_copy_to(ptr(), to);
}
inline void Matrix::update(const Matrix &from) {
  SafeCall(igraph_matrix_update(ptr(), from.ptr()));
}
inline void Matrix::swap(Matrix &m2) {
  SafeCall(igraph_matrix_swap(ptr(), m2.ptr()));
}

/* Accessing elements */
inline double &Matrix::operator()(long int i, long int j) noexcept {
  return MATRIX(*ptr(), i, j);
}

/* Operations on rows and columns */
inline Vector Matrix::get_row(long int index) const {
  Vector vector(ncol());
  SafeCall(igraph_matrix_get_row(ptr(), vector.ptr(), index));
  return vector;
}
inline Vector Matrix::get_col(long int index) const {
  Vector vector(nrow());
  SafeCall(igraph_matrix_get_col(ptr(), vector.ptr(), index));
  return vector;
}
inline void Matrix::set_row(const Vector &vector, long int index) {
  SafeCall(igraph_matrix_set_row(ptr(), vector.ptr(), index));
}
inline void Matrix::set_col(const Vector &vector, long int index) {
  SafeCall(igraph_matrix_set_col(ptr(), vector.ptr(), index));
}
inline void Matrix::swap_rows(long int i, long j) {
  SafeCall(igraph_matrix_swap_rows(ptr(), i, j));
}
inline void Matrix::swap_cols(long int i, long j) {
  SafeCall(igraph_matrix_swap_cols(ptr(), i, j));
}
inline Matrix Matrix::select_rows(const Vector &rows) {
  Matrix matrix(rows.size(), ncol());
  SafeCall(igraph_matrix_select_rows(ptr(), matrix.ptr(), rows.ptr()));
  return matrix;
}
inline Matrix Matrix::select_cols(const Vector &cols) {
  Matrix matrix(nrow(), cols.size());
  SafeCall(igraph_matrix_select_cols(ptr(), matrix.ptr(), cols.ptr()));
  return matrix;
}
inline Matrix Matrix::select_rows_cols(const Vector &rows, const Vector &cols) {
  Matrix matrix(rows.size(), cols.size());
  SafeCall(igraph_matrix_select_rows_cols(ptr(), matrix.ptr(), rows.ptr(),
                                          cols.ptr()));
  return matrix;
}

/* Matrix operations */
inline Matrix &Matrix::operator+=(double scalar) noexcept {
  igraph_matrix_add_constant(ptr(), scalar);
  return *this;
}
inline Matrix &Matrix::operator-=(double scalar) noexcept {
  return this->operator+=(-scalar);
}
inline Matrix &Matrix::operator*=(double scalar) noexcept {
  igraph_matrix_scale(ptr(), scalar);
  return *this;
}
inline Matrix &Matrix::operator/=(double scalar) noexcept {
  return this->operator*=(1.0 / scalar);
}
inline Matrix &Matrix::operator+=(const Matrix &matrix) {
  SafeCall(igraph_matrix_add(ptr(), matrix.ptr()));
  return *this;
}
inline Matrix &Matrix::operator-=(const Matrix &matrix) {
  SafeCall(igraph_matrix_sub(ptr(), matrix.ptr()));
  return *this;
}
inline Matrix &Matrix::operator*=(const Matrix &matrix) {
  SafeCall(igraph_matrix_mul_elements(ptr(), matrix.ptr()));
  return *this;
}
inline Matrix &Matrix::operator/=(const Matrix &matrix) {
  SafeCall(igraph_matrix_div_elements(ptr(), matrix.ptr()));
  return *this;
}
inline double Matrix::sum() const noexcept { return igraph_matrix_sum(ptr()); }
inline double Matrix::prod() const noexcept {
  return igraph_matrix_prod(ptr());
}
inline Vector Matrix::rowsums() const {
  Vector vector(ncol());
  SafeCall(igraph_matrix_rowsum(ptr(), vector.ptr()));
  return vector;
}
inline Vector Matrix::colsums() const {
  Vector vector(nrow());
  SafeCall(igraph_matrix_colsum(ptr(), vector.ptr()));
  return vector;
}
inline void Matrix::transpose() { SafeCall(igraph_matrix_transpose(ptr())); }

/* Matrix comparisons */
inline bool Matrix::operator==(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_e(ptr(), matrix.ptr());
}
inline bool Matrix::operator<(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_l(ptr(), matrix.ptr());
}
inline bool Matrix::operator>(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_g(ptr(), matrix.ptr());
}
inline bool Matrix::operator<=(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_le(ptr(), matrix.ptr());
}
inline bool Matrix::operator>=(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_ge(ptr(), matrix.ptr());
}

/* Combining matrices */
inline void Matrix::rbind(const Matrix &matrix) {
  SafeCall(igraph_matrix_rbind(ptr(), matrix.ptr()));
}
inline void Matrix::cbind(const Matrix &matrix) {
  SafeCall(igraph_matrix_cbind(ptr(), matrix.ptr()));
}

/* Finding minimum and maximum */
inline double Matrix::min() const noexcept { return igraph_matrix_min(ptr()); }
inline double Matrix::max() const noexcept { return igraph_matrix_max(ptr()); }
inline void Matrix::which_min(long int *i, long int *j) const {
  SafeCall(igraph_matrix_which_min(ptr(), i, j));
}
inline void Matrix::which_max(long int *i, long int *j) const {
  SafeCall(igraph_matrix_which_max(ptr(), i, j));
}
inline void Matrix::minmax(double *min, double *max) const {
  SafeCall(igraph_matrix_minmax(ptr(), min, max));
}
inline void Matrix::which_minmax(long int *imin, long int *jmin, long int *imax,
                                 long int *jmax) const {
  SafeCall(igraph_matrix_which_minmax(ptr(), imin, jmin, imax, jmax));
}

/* Matrix properties */
inline bool Matrix::empty() const noexcept {
  return igraph_matrix_empty(ptr());
}
inline bool Matrix::isnull() const noexcept {
  return igraph_matrix_isnull(ptr());
}
inline long int Matrix::size() const noexcept {
  return igraph_matrix_size(ptr());
}
inline long int Matrix::capacity() const noexcept {
  return igraph_matrix_capacity(ptr());
}
inline long int Matrix::nrow() const noexcept {
  return igraph_matrix_nrow(ptr());
}
inline long int Matrix::ncol() const noexcept {
  return igraph_matrix_ncol(ptr());
}
inline bool Matrix::is_symmetric() const noexcept {
  return igraph_matrix_is_symmetric(ptr());
}
inline double Matrix::maxdifference(const Matrix &matrix) const noexcept {
  return igraph_matrix_maxdifference(ptr(), matrix.ptr());
}

/* Searching for elements */
inline bool Matrix::contains(double value) const noexcept {
  return igraph_matrix_contains(ptr(), value);
}
inline bool Matrix::search(double value, long int *row, long int *col,
                           long int *pos, double from_pos) const noexcept {
  return igraph_matrix_search(ptr(), from_pos, value, pos, row, col);
}

/* Resize operations */
inline void Matrix::resize(long int nrow, long int ncol) {
  SafeCall(igraph_matrix_resize(ptr(), nrow, ncol));
}
inline void Matrix::resize_min() { SafeCall(igraph_matrix_resize_min(ptr())); }
inline void Matrix::add_rows(long int n) {
  SafeCall(igraph_matrix_add_rows(ptr(), n));
}
inline void Matrix::add_cols(long int n) {
  SafeCall(igraph_matrix_add_cols(ptr(), n));
}
inline void Matrix::remove_row(long int row) {
  SafeCall(igraph_matrix_remove_row(ptr(), row));
}
inline void Matrix::remove_col(long int col) {
  SafeCall(igraph_matrix_remove_col(ptr(), col));
}

} // namespace igraph

#endif // IGRAPHPP_MATRIX_IMPL_HPP_
