
find_package(Catch REQUIRED)
include_directories(${CATCH_INCLUDE_DIR})

add_executable(test_suite
	test_suite.cpp
	VectorView.cpp)

target_link_libraries(test_suite ${IGRAPH_LIBRARIES})
