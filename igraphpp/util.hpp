#ifndef IGRAPHPP_UTIL_HPP_
#define IGRAPHPP_UTIL_HPP_

#ifndef IGRAPHPP_IGRAPH_HPP_
#error "You must include igraph.hpp first"
#endif

#include <iterator>
#include <type_traits>

namespace igraph {

namespace util {

template <typename T, typename = void>
struct is_iterator {
  static constexpr bool value = false;
};

template <typename T>
struct is_iterator<
    T, typename std::enable_if<!std::is_same<
           typename std::iterator_traits<T>::value_type, void>::value>::type> {
  static constexpr bool value = true;
};

constexpr bool all_args() { return true; }

template <typename... Tail>
constexpr bool all_args(bool head, Tail... tail) {
  return head && all_args(tail...);
}

}  // namespace util

}  // namespace igraph

#endif  // IGRAPHPP_UTIL_HPP_
