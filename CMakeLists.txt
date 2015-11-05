
cmake_minimum_required(VERSION 2.8.9)
project(IGraphWrapper)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMakeModules)

find_package(IGRAPH REQUIRED)

include_directories(${IGRAPH_INCLUDES})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra")

add_subdirectory(examples)
add_subdirectory(test)
