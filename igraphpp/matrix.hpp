#ifndef IGRAPHPP_MATRIX_HPP_
#define IGRAPHPP_MATRIX_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <initializer_list>

#include <igraph.h>

#include "./mapper.hpp"
#include "./util.hpp"
#include "./vectorptr.hpp"
#include "./vector.hpp"

namespace igraph {

class Vector;

class Matrix {
  template <typename>
  friend struct TypeMapper;

 public:
  /* Constructors and destructors */
  ~Matrix();
  Matrix(long int nrow = 0, long int ncol = 0);
  Matrix(std::initializer_list<double> elements, long int nrow);
  Matrix(const VectorView &elements, long int nrow);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  Matrix(Iterator begin, Iterator end, long int nrow, long int ncol);
  Matrix(const Matrix &matrix);
  Matrix(Matrix &&matrix);
  Matrix &operator=(const Matrix &matrix)&;
  Matrix &operator=(Matrix &&matrix)&;

  /* Initializing elements */
  void null() noexcept;
  void fill(double value) noexcept;

  /* Copying matrices */
  void copy_to(double *to) noexcept;  // Column major
  void update(const Matrix &from);
  void swap(Matrix &m2);

  /* Accessing elements */
  double &operator()(long int i, long int j) noexcept;
  double &at(long int i, long int j) noexcept;

  /* Operations on rows and columns */
  Vector get_row(long int index) const;
  Vector get_col(long int index) const;
  void set_row(long int index, const Vector &vector);
  void set_col(long int index, const Vector &vector);
  void swap_rows(long int i, long j);
  void swap_cols(long int i, long j);
  Matrix select_rows(const Vector &rows);
  Matrix select_cols(const Vector &cols);
  Matrix select_rows_cols(const Vector &rows, const Vector &cols);

  /* Matrix operations */
  Matrix &operator+=(double scalar) noexcept;
  Matrix &operator-=(double scalar) noexcept;
  Matrix &operator*=(double scalar) noexcept;
  Matrix &operator/=(double scalar) noexcept;
  Matrix &operator+=(const Matrix &matrix);
  Matrix &operator-=(const Matrix &matrix);
  Matrix &operator*=(const Matrix &matrix);
  Matrix &operator/=(const Matrix &matrix);
  double sum() const noexcept;
  double prod() const noexcept;
  Vector rowsums() const;
  Vector colsums() const;
  void transpose();

  /* Matrix comparisons */
  bool operator==(const Matrix &matrix) const noexcept;
  bool operator<(const Matrix &matrix) const noexcept;
  bool operator>(const Matrix &matrix) const noexcept;
  bool operator<=(const Matrix &matrix) const noexcept;
  bool operator>=(const Matrix &matrix) const noexcept;

  /* Combining matrices */
  void rbind(const Matrix &matrix);
  void cbind(const Matrix &matrix);

  /* Finding minimum and maximum */
  double min() const noexcept;
  double max() const noexcept;
  void which_min(long int *i, long int *j) const;
  void which_max(long int *i, long int *j) const;
  void minmax(double *min, double *max) const;
  void which_minmax(long int *imin, long int *jmin, long int *imax,
                    long int *jmax) const;

  /* Matrix properties */
  bool empty() const noexcept;
  bool isnull() const noexcept;
  long int size() const noexcept;
  long int capacity() const noexcept;
  long int nrow() const noexcept;
  long int ncol() const noexcept;
  bool is_symmetric() const noexcept;
  double maxdifference(const Matrix &matrix) const noexcept;

  /* Searching for elements */
  bool contains(double value) const noexcept;
  bool search(double value, long int *row, long int *col, long int *pos = NULL,
              double from_pos = 0) const noexcept;

  /* Resize operations */
  void resize(long int nrow, long int ncol);
  void resize_min();
  void add_rows(long int n);
  void add_cols(long int n);
  void remove_row(long int row);
  void remove_col(long int col);

  igraph_matrix_t *ptr() { return &matrix_; }
  const igraph_matrix_t *ptr() const { return &matrix_; }

 protected:
  Matrix(const igraph_matrix_t &matrix, bool owner = true);

 private:
  void disown() { owner_ = false; }
  bool owner() const { return owner_; }

  igraph_matrix_t matrix_;
  bool owner_ = true;
};

template <>
struct TypeMapper<Matrix> {
  typedef igraph_matrix_t type;

  static Matrix Build(const type &vector) { return Matrix(vector, false); }
};

template <typename Iterator, typename>
Matrix::Matrix(Iterator begin, Iterator end, long int nrow, long int ncol)
    : Matrix(nrow, ncol) {
  long int i = 0, j = 0;
  for (; begin != end; ++begin) {
    this->at(i, j) = *begin;
    ++i;
    if (i == nrow) {
      i = 0;
      ++j;
    }
  }
}

}  // namespace igraph

#endif  // IGRAPHPP_MATRIX_HPP_
