#ifndef __IGRAPHPP_EXCEPTION_HPP_
#define __IGRAPHPP_EXCEPTION_HPP_

#include <stdexcept>

#include <igraph.h>

namespace igraph {

class Exception : public std::exception {
public:
  explicit Exception(int code) : std::exception(), errno_(code){};

  virtual ~Exception() throw() {}

  virtual const char *what() const throw() {
    return igraph_strerror(static_cast<int>(errno_));
  }

private:
  int errno_;
};

static inline int SafeCall(int ret) {
  if (ret != IGRAPH_SUCCESS)
    throw Exception(ret);
  return ret;
}

} // namespace igraph

#endif // __IGRAPHPP_EXCEPTION_HPP_
