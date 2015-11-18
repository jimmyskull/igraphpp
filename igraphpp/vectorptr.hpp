#ifndef IGRAPH_VECTORPTR_HPP_
#define IGRAPH_VECTORPTR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph/igraph.h>

#include "./exception.hpp"
#include "./mapper.hpp"

namespace igraph {

template <typename Type>
class VectorPtr {
 private:
  typedef TypeMapper<Type> mapper;
  typedef Type wrapper_t;
  typedef typename mapper::type wrapped_t;

 public:
  /* Constructors and Destructors */
  ~VectorPtr() { igraph_vector_ptr_destroy(ptr()); }
  VectorPtr(long int size = 0) {
    SafeCall(igraph_vector_ptr_init(ptr(), size));
  }
  VectorPtr(const VectorPtr<Type> &other) {
    SafeCall(igraph_vector_ptr_copy(ptr(), other.ptr()));
  }

  void free_all() noexcept { igraph_vector_ptr_free_all(ptr()); }

  long int size() const noexcept { return igraph_vector_ptr_size(ptr()); }

  void clear() noexcept { igraph_vector_ptr_clear(ptr()); }

  void push_back(wrapper_t &element) {
    void *element_pointer = reinterpret_cast<void *>(element.ptr());
    SafeCall(igraph_vector_ptr_push_back(ptr(), element_pointer));
  }

  wrapper_t operator[](long int pos) const {
    void *element = VECTOR(*ptr())[pos];
    return mapper::Build(*reinterpret_cast<wrapped_t *>(element));
  }

  void resize(long int newsize) {
    SafeCall(igraph_vector_ptr_resize(ptr(), newsize));
  }

  igraph_vector_ptr_t *ptr() { return &vector_; }
  const igraph_vector_ptr_t *ptr() const { return &vector_; }

 private:
  igraph_vector_ptr_t vector_;
};

}  // namespace igraph

#endif  // IGRAPH_VECTORPTR_HPP_
