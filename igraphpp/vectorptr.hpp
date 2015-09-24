#ifndef IGRAPHPP_VECTORPTR_HPP_
#define IGRAPHPP_VECTORPTR_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <igraph.h>

namespace igraph {

template <typename PtrType = void *, typename PtrWrapper = void *>
class VectorPtr {
public:
  /* Constructors and Destructors */
  ~VectorPtr() { igraph_vector_ptr_destroy(ptr()); }
  explicit VectorPtr(long int size) {
    SafeCall(igraph_vector_ptr_init(ptr(), size));
  }
  VectorPtr(const VectorPtr &other) {
    SafeCall(igraph_vector_ptr_copy(ptr(), other.ptr()));
  }

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
  }

  igraph_vector_ptr_t *ptr() { return &vector_; }
  const igraph_vector_ptr_t *ptr() const { return &vector_; }

private:
  igraph_vector_ptr_t vector_;
};

} // namespace igraph

#endif // IGRAPHPP_VECTORPTR_HPP_
