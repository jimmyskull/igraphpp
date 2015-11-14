#ifndef IGRAPHPP_VECTOR_HPP_
#define IGRAPHPP_VECTOR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <initializer_list>
#include <type_traits>

#include <igraph.h>

#include "./exception.hpp"
#include "./util.hpp"

namespace igraph {

class Vector;

class VectorView {
 public:
  /* Constructors and Destructors */
  ~VectorView();
  VectorView(const double *data, long int length);
  VectorView &operator=(VectorView &other) = delete;
  VectorView &operator=(VectorView &&other) = delete;

  /* Initializing elements */
  void null() noexcept;
  void fill(double value) noexcept;

  /* Accessing elements */
  double &operator[](long int i) const noexcept;
  double &at(long int i) const noexcept;
  double &head() const noexcept;
  double &tail() const noexcept;

  /* Copying vectors */
  void copy_to(double *output) noexcept;

  /* Exchanging elements */
  void swap(long int i, long int j) noexcept;
  void swap(VectorView &b);
  void shuffle() noexcept;
  void reverse() noexcept;

  /* Vector operations */
  VectorView &operator+=(double scalar) noexcept;
  VectorView &operator-=(double scalar) noexcept;
  VectorView &operator*=(double scalar) noexcept;
  VectorView &operator/=(double scalar) noexcept;
  VectorView &operator+=(const VectorView &b);
  VectorView &operator-=(const VectorView &b);
  VectorView &operator*=(const VectorView &b);
  VectorView &operator/=(const VectorView &b);

  /* Vector comparisons */
  bool operator==(const VectorView &b) const noexcept;
  bool operator!=(const VectorView &b) const noexcept;
  bool operator<(const VectorView &b) const noexcept;
  bool operator>(const VectorView &b) const noexcept;
  bool operator<=(const VectorView &b) const noexcept;
  bool operator>=(const VectorView &b) const noexcept;

  /* Finding minimum and maximum */
  double min() const noexcept;
  double max() const noexcept;
  long int which_min() const noexcept;
  long int which_max() const noexcept;
  void minmax(double *min, double *max) const;
  void which_minmax(long int *min, long int *max) const;
  /* Vector properties */
  bool empty() const noexcept;
  long int size() const noexcept;
  long int capacity() const noexcept;
  double sum() const noexcept;
  double prod() const noexcept;
  bool isininterval(double low, double high) const noexcept;
  double maxdifference(const VectorView &b) const noexcept;
  /* Searching for elements */
  bool contains(double value) const noexcept;
  long int search(double value, long int from = 0) const noexcept;
  long int binsearch(double value) const noexcept;

  /* Sorting */
  void sort() noexcept;

  /* Internal use */
  const igraph_vector_t *ptr() const { return is_none() ? NULL : &vector_; }
  igraph_vector_t *ptr() { return is_none() ? NULL : &vector_; }

  class iterator : public std::iterator<std::forward_iterator_tag, double> {
   public:
    iterator(igraph_vector_t *vec, long int pos) : vec_(vec), pos_(pos) {}
    iterator(const iterator &itr) : vec_(itr.vec_), pos_(itr.pos_) {}
    iterator &operator++() {
      ++pos_;
      return *this;
    }
    bool operator==(const iterator &rhs) const { return pos_ == rhs.pos_; }
    bool operator!=(const iterator &rhs) const { return pos_ != rhs.pos_; }
    double &operator*() { return VECTOR(*vec_)[pos_]; }

   private:
    igraph_vector_t *vec_;
    long int pos_;
  };

  class const_iterator
      : public std::iterator<std::forward_iterator_tag, double> {
   public:
    const_iterator(const igraph_vector_t *vec, long int pos)
        : vec_(vec), pos_(pos) {}
    const_iterator(const const_iterator &itr)
        : vec_(itr.vec_), pos_(itr.pos_) {}
    const_iterator &operator++() {
      ++pos_;
      return *this;
    }
    bool operator==(const const_iterator &rhs) { return pos_ == rhs.pos_; }
    bool operator!=(const const_iterator &rhs) { return pos_ != rhs.pos_; }
    double operator*() const { return VECTOR(*vec_)[pos_]; }

   private:
    const igraph_vector_t *vec_;
    long int pos_;
  };

  iterator begin() { return iterator(ptr(), 0); }
  iterator end() { return iterator(ptr(), size()); }

  const_iterator cbegin() const { return const_iterator(ptr(), 0); }
  const_iterator cend() const { return const_iterator(ptr(), size()); }

  static VectorView None() {
    static VectorView instance = VectorView(true);
    return instance;
  }
  bool is_none() const { return none_vector_; }

 protected:
  VectorView() = default;
  VectorView(bool none_vector) : none_vector_(none_vector) {}

 private:
  VectorView(const VectorView &) = default;

  igraph_vector_t vector_;
  bool none_vector_ = false;
};

class Vector : public VectorView {
 public:
  /* Constructors and Destructors */
  ~Vector();
  explicit Vector(long int size = 0);
  Vector(const double *data, long int length);
  Vector(double from, double to);
  Vector(std::initializer_list<double> list);
  template <typename Iterator, typename = typename std::enable_if<
                                   util::is_iterator<Iterator>::value>::type>
  Vector(Iterator begin, Iterator end);
  Vector(const Vector &other);
  Vector(const VectorView &other);
  Vector(Vector &&other);
  Vector &operator=(Vector &&other);

  /* Vector operations */
  Vector operator+(double scalar) const noexcept;
  Vector operator-(double scalar) const noexcept;
  Vector operator*(double scalar) const noexcept;
  Vector operator/(double scalar) const noexcept;
  Vector operator+(const VectorView &b) const;
  Vector operator-(const VectorView &b) const;
  Vector operator*(const VectorView &b) const;
  Vector operator/(const VectorView &b) const;

  /* Copying vectors */
  Vector &operator=(const VectorView &other);
  Vector &operator=(const Vector &other);
  void append(const VectorView &other);
  void update(const VectorView &other) noexcept;

  /* Resizing operations */
  void clear() noexcept;
  void reserve(long int newsize);
  void resize(long int newsize);
  void resize_min();
  void push_back(double value);
  double pop_back() noexcept;
  void insert(long int pos, double value);
  void remove(long int pos) noexcept;
  void remove(long int from, long int to) noexcept;

  /* Set operations on sorted vectors */
  Vector intersect_sorted(const VectorView &other) const;
  Vector difference_sorted(const VectorView &other) const;

  static Vector Repeat(double value, long int times);

 private:
  explicit Vector(const igraph_vector_t &vector);

  void disown() { VECTOR(*ptr()) = NULL; }
  bool owner() const { return VECTOR(*ptr()) != NULL; }
};

template <typename Iterator, typename>
Vector::Vector(Iterator begin, Iterator end)
    : VectorView() {
  SafeCall(igraph_vector_init(ptr(), 0));
  for (; begin != end; ++begin) {
    push_back(static_cast<double>(*begin));
  }
}

}  // namespace igraph

#endif  // IGRAPHPP_VECTOR_HPP_
