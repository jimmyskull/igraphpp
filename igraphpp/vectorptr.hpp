<<<<<<< HEAD
#ifndef IGRAPHPP_VECTORPTR_HPP_
#define IGRAPHPP_VECTORPTR_HPP_
=======
#ifndef IGRAPH_VECTORPTR_HPP_
#define IGRAPH_VECTORPTR_HPP_
>>>>>>> 2ac277531cb77475e8e38bdd4438ac500a5d5628

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

<<<<<<< HEAD
namespace igraph {

template <typename PtrType = void *, typename PtrWrapper = void *>
class VectorPtr {
public:
  /* Constructors and Destructors */
  ~VectorPtr() { igraph_vector_ptr_destroy(ptr()); }
  explicit VectorPtr(long int size) {
=======
#include "./exception.hpp"

namespace igraph {

template <typename PtrType = void *> class VectorPtr {
public:
  /* Constructors and Destructors */
  ~VectorPtr() { igraph_vector_ptr_destroy(ptr()); }
  VectorPtr(long int size = 0) {
>>>>>>> 2ac277531cb77475e8e38bdd4438ac500a5d5628
    SafeCall(igraph_vector_ptr_init(ptr(), size));
  }
  VectorPtr(const VectorPtr &other) {
    SafeCall(igraph_vector_ptr_copy(ptr(), other.ptr()));
  }

<<<<<<< HEAD
  void free_all() { igraph_vector_ptr_free_all(ptr()); }
  long int size() const noexcept { return igraph_vector_ptr_size(ptr()); }
  void clear() noexcept { igraph_vector_ptr_clear(ptr()); }
  void push_back(PtrType element) {
    SafeCall(igraph_vector_ptr_push_back(ptr(),
          reinterpret_cast<void *>(element)));
  }

  PtrWrapper operator[](long int pos) const noexcept {
    PtrType element_ptr = reinterpret_cast<PtrType>(VECTOR(*ptr())[pos]);
    return PtrWrapper(element_ptr);
=======
  void free_all() noexcept { igraph_vector_ptr_free_all(ptr()); }

  long int size() const noexcept { return igraph_vector_ptr_size(ptr()); }

  void clear() noexcept { igraph_vector_ptr_clear(ptr()); }

  void push_back(PtrType element) {
    void *element_pointer = reinterpret_cast<void *>(element);
    SafeCall(igraph_vector_ptr_push_back(ptr(), element_pointer));
  }

  PtrType operator[](long int pos) const {
    void *element = VECTOR(vector())[pos];
    return reinterpret_cast<PtrType>(element);
  }

  void resize(long int newsize) {
    SafeCall(igraph_vector_ptr_resize(ptr(), newsize));
>>>>>>> 2ac277531cb77475e8e38bdd4438ac500a5d5628
  }

  igraph_vector_ptr_t *ptr() { return &vector_; }
  const igraph_vector_ptr_t *ptr() const { return &vector_; }

private:
<<<<<<< HEAD
=======
  const igraph_vector_ptr_t &vector() const { return vector_; }

>>>>>>> 2ac277531cb77475e8e38bdd4438ac500a5d5628
  igraph_vector_ptr_t vector_;
};

} // namespace igraph

<<<<<<< HEAD
#endif // IGRAPHPP_VECTORPTR_HPP_
=======
#endif // IGRAPH_VECTORPTR_HPP_
>>>>>>> 2ac277531cb77475e8e38bdd4438ac500a5d5628
