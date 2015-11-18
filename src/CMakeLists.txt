
add_library(igraphpp STATIC graph.cpp matrix.cpp edge_iterator.cpp
  edge_selector.cpp vertex_iterator.cpp vertex_selector.cpp vector.cpp)
target_link_libraries(igraphpp ${IGRAPH_LIBRARIES})
