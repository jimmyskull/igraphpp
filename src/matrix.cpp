#include "../igraphpp/igraph.hpp"

#include <igraph.h>

namespace igraph {

/* Constructors and destructors */
Matrix::~Matrix() {
  if (owner()) igraph_matrix_destroy(ptr());
}
Matrix::Matrix(long int nrow, long int ncol) {
  SafeCall(igraph_matrix_init(ptr(), nrow, ncol));
}
Matrix::Matrix(std::initializer_list<double> elements, long int nrow)
    : Matrix(elements.begin(), elements.end(), nrow, elements.size() / nrow) {}
Matrix::Matrix(const VectorView &elements, long int nrow)
    : Matrix(elements.cbegin(), elements.cend(), nrow, elements.size() / nrow) {
}
Matrix::Matrix(const Matrix &matrix) {
  SafeCall(igraph_matrix_copy(ptr(), matrix.ptr()));
}
Matrix::Matrix(Matrix &&matrix) {
  matrix_ = *matrix.ptr();
  matrix.disown();
}
Matrix &Matrix::operator=(const Matrix &matrix) & {
  this->~Matrix();
  SafeCall(igraph_matrix_copy(ptr(), matrix.ptr()));
  return *this;
}
Matrix &Matrix::operator=(Matrix &&matrix) & {
  this->~Matrix();
  *ptr() = *matrix.ptr();
  matrix.disown();
  return *this;
}

/* Initializing elements */
void Matrix::null() noexcept { igraph_matrix_null(ptr()); }
void Matrix::fill(double value) noexcept { igraph_matrix_fill(ptr(), value); }

/* Copying matrices */
void Matrix::copy_to(double *to) noexcept { igraph_matrix_copy_to(ptr(), to); }
void Matrix::update(const Matrix &from) {
  SafeCall(igraph_matrix_update(ptr(), from.ptr()));
}
void Matrix::swap(Matrix &m2) { SafeCall(igraph_matrix_swap(ptr(), m2.ptr())); }

/* Accessing elements */
double &Matrix::operator()(long int i, long int j) noexcept {
  return MATRIX(*ptr(), i, j);
}
double &Matrix::at(long int i, long int j) noexcept {
  return MATRIX(*ptr(), i, j);
}

/* Operations on rows and columns */
Vector Matrix::get_row(long int index) const {
  Vector vector(ncol());
  SafeCall(igraph_matrix_get_row(ptr(), vector.ptr(), index));
  return vector;
}
Vector Matrix::get_col(long int index) const {
  Vector vector(nrow());
  SafeCall(igraph_matrix_get_col(ptr(), vector.ptr(), index));
  return vector;
}
void Matrix::set_row(long int index, const Vector &vector) {
  SafeCall(igraph_matrix_set_row(ptr(), vector.ptr(), index));
}
void Matrix::set_col(long int index, const Vector &vector) {
  SafeCall(igraph_matrix_set_col(ptr(), vector.ptr(), index));
}
void Matrix::swap_rows(long int i, long j) {
  SafeCall(igraph_matrix_swap_rows(ptr(), i, j));
}
void Matrix::swap_cols(long int i, long j) {
  SafeCall(igraph_matrix_swap_cols(ptr(), i, j));
}
Matrix Matrix::select_rows(const Vector &rows) {
  Matrix matrix(rows.size(), ncol());
  SafeCall(igraph_matrix_select_rows(ptr(), matrix.ptr(), rows.ptr()));
  return matrix;
}
Matrix Matrix::select_cols(const Vector &cols) {
  Matrix matrix(nrow(), cols.size());
  SafeCall(igraph_matrix_select_cols(ptr(), matrix.ptr(), cols.ptr()));
  return matrix;
}
Matrix Matrix::select_rows_cols(const Vector &rows, const Vector &cols) {
  Matrix matrix(rows.size(), cols.size());
  SafeCall(igraph_matrix_select_rows_cols(ptr(), matrix.ptr(), rows.ptr(),
                                          cols.ptr()));
  return matrix;
}

/* Matrix operations */
Matrix &Matrix::operator+=(double scalar) noexcept {
  igraph_matrix_add_constant(ptr(), scalar);
  return *this;
}
Matrix &Matrix::operator-=(double scalar) noexcept {
  return this->operator+=(-scalar);
}
Matrix &Matrix::operator*=(double scalar) noexcept {
  igraph_matrix_scale(ptr(), scalar);
  return *this;
}
Matrix &Matrix::operator/=(double scalar) noexcept {
  return this->operator*=(1.0 / scalar);
}
Matrix &Matrix::operator+=(const Matrix &matrix) {
  SafeCall(igraph_matrix_add(ptr(), matrix.ptr()));
  return *this;
}
Matrix &Matrix::operator-=(const Matrix &matrix) {
  SafeCall(igraph_matrix_sub(ptr(), matrix.ptr()));
  return *this;
}
Matrix &Matrix::operator*=(const Matrix &matrix) {
  SafeCall(igraph_matrix_mul_elements(ptr(), matrix.ptr()));
  return *this;
}
Matrix &Matrix::operator/=(const Matrix &matrix) {
  SafeCall(igraph_matrix_div_elements(ptr(), matrix.ptr()));
  return *this;
}
double Matrix::sum() const noexcept { return igraph_matrix_sum(ptr()); }
double Matrix::prod() const noexcept { return igraph_matrix_prod(ptr()); }
Vector Matrix::rowsums() const {
  Vector vector(ncol());
  SafeCall(igraph_matrix_rowsum(ptr(), vector.ptr()));
  return vector;
}
Vector Matrix::colsums() const {
  Vector vector(nrow());
  SafeCall(igraph_matrix_colsum(ptr(), vector.ptr()));
  return vector;
}
void Matrix::transpose() { SafeCall(igraph_matrix_transpose(ptr())); }

/* Matrix comparisons */
bool Matrix::operator==(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_e(ptr(), matrix.ptr());
}
bool Matrix::operator<(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_l(ptr(), matrix.ptr());
}
bool Matrix::operator>(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_g(ptr(), matrix.ptr());
}
bool Matrix::operator<=(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_le(ptr(), matrix.ptr());
}
bool Matrix::operator>=(const Matrix &matrix) const noexcept {
  return igraph_matrix_all_ge(ptr(), matrix.ptr());
}

/* Combining matrices */
void Matrix::rbind(const Matrix &matrix) {
  SafeCall(igraph_matrix_rbind(ptr(), matrix.ptr()));
}
void Matrix::cbind(const Matrix &matrix) {
  SafeCall(igraph_matrix_cbind(ptr(), matrix.ptr()));
}

/* Finding minimum and maximum */
double Matrix::min() const noexcept { return igraph_matrix_min(ptr()); }
double Matrix::max() const noexcept { return igraph_matrix_max(ptr()); }
void Matrix::which_min(long int *i, long int *j) const {
  SafeCall(igraph_matrix_which_min(ptr(), i, j));
}
void Matrix::which_max(long int *i, long int *j) const {
  SafeCall(igraph_matrix_which_max(ptr(), i, j));
}
void Matrix::minmax(double *min, double *max) const {
  SafeCall(igraph_matrix_minmax(ptr(), min, max));
}
void Matrix::which_minmax(long int *imin, long int *jmin, long int *imax,
                          long int *jmax) const {
  SafeCall(igraph_matrix_which_minmax(ptr(), imin, jmin, imax, jmax));
}

/* Matrix properties */
bool Matrix::empty() const noexcept { return igraph_matrix_empty(ptr()); }
bool Matrix::isnull() const noexcept { return igraph_matrix_isnull(ptr()); }
long int Matrix::size() const noexcept { return igraph_matrix_size(ptr()); }
long int Matrix::capacity() const noexcept {
  return igraph_matrix_capacity(ptr());
}
long int Matrix::nrow() const noexcept { return igraph_matrix_nrow(ptr()); }
long int Matrix::ncol() const noexcept { return igraph_matrix_ncol(ptr()); }
bool Matrix::is_symmetric() const noexcept {
  return igraph_matrix_is_symmetric(ptr());
}
double Matrix::maxdifference(const Matrix &matrix) const noexcept {
  return igraph_matrix_maxdifference(ptr(), matrix.ptr());
}

/* Searching for elements */
bool Matrix::contains(double value) const noexcept {
  return igraph_matrix_contains(ptr(), value);
}
bool Matrix::search(double value, long int *row, long int *col, long int *pos,
                    double from_pos) const noexcept {
  return igraph_matrix_search(ptr(), from_pos, value, pos, row, col);
}

/* Resize operations */
void Matrix::resize(long int nrow, long int ncol) {
  SafeCall(igraph_matrix_resize(ptr(), nrow, ncol));
}
void Matrix::resize_min() { SafeCall(igraph_matrix_resize_min(ptr())); }
void Matrix::add_rows(long int n) {
  SafeCall(igraph_matrix_add_rows(ptr(), n));
}
void Matrix::add_cols(long int n) {
  SafeCall(igraph_matrix_add_cols(ptr(), n));
}
void Matrix::remove_row(long int row) {
  SafeCall(igraph_matrix_remove_row(ptr(), row));
}
void Matrix::remove_col(long int col) {
  SafeCall(igraph_matrix_remove_col(ptr(), col));
}

Matrix::Matrix(const igraph_matrix_t &matrix, bool owner)
    : matrix_(matrix), owner_(owner) {}

}  // namespace igraph
