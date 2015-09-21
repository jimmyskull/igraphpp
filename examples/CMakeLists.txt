
add_executable(helloworld hello_world.cpp)
target_link_libraries(helloworld ${IGRAPH_LIBRARIES})

add_executable(vector vector.cpp)
target_link_libraries(vector ${IGRAPH_LIBRARIES})

add_executable(erdosrenyi erdos_renyi.cpp)
target_link_libraries(erdosrenyi ${IGRAPH_LIBRARIES})

