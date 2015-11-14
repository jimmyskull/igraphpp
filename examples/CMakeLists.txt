
add_executable(helloworld hello_world.cpp)
target_link_libraries(helloworld ${IGRAPH_LIBRARIES} igraphpp)

add_executable(vector vector.cpp)
target_link_libraries(vector ${IGRAPH_LIBRARIES} igraphpp)

add_executable(lesson1 lesson1.cpp)
target_link_libraries(lesson1 ${IGRAPH_LIBRARIES} igraphpp)

add_executable(lesson2 lesson2.cpp)
target_link_libraries(lesson2 ${IGRAPH_LIBRARIES} igraphpp)

add_executable(lesson3 lesson3.cpp)
target_link_libraries(lesson3 ${IGRAPH_LIBRARIES} igraphpp)
