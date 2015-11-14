cmake_minimum_required(VERSION 2.8.9)
project(IGraphWrapper C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMakeModules)

find_package(IGRAPH REQUIRED)

include_directories(${IGRAPH_INCLUDES})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra")

Option(BuildExamples "BuildExamples" OFF)
Option(BuildTest "BuildTest" OFF)

add_subdirectory(src)

if (BuildExamples)
  add_subdirectory(examples)
endif (BuildExamples)

if (BuildTest)
  add_subdirectory(test)
endif (BuildTest)

