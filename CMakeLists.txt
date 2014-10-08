cmake_minimum_required(VERSION 2.8)

find_package(PythonLibs)
message(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
set(numpy_include "/home/js850/.local/lib/python2.7/site-packages/numpy/core/include")
include_directories(${numpy_include})

enable_language(CXX)
SET(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -std=c++0x")

#cmake_policy(SET CMP0015 NEW)

# set the mcpele include directory
set(mcpele_include ${CMAKE_SOURCE_DIR}/source)
include_directories(${mcpele_include})
message("mcpele include directory ${mcpele_include}")

# set the pele include directory
# the pele include directory should be specified from python
set(pele_include ${CMAKE_SOURCE_DIR}/../pele/source CACHE STRING "the pele c++ source directory")
if (EXISTS ${pele_include}/pele/array.h)
  message("pele include directory: ${pele_include}")
else()
  message(FATAL_ERROR "pele include directory is not correct: ${pele_include} : use ccmake to set it")
endif()
include_directories(${pele_include})
FILE(GLOB pele_sources ${pele_include}/*.cpp)

# build the mcpele library
FILE(GLOB mcpele_sources ${mcpele_include}/*.cpp)
add_library(mcpele_lib SHARED ${pele_sources} ${mcpele_sources})

set(action_source ${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_action_cpp.cxx)
#add_library(_action_cpp SHARED ${action_source} ${pele_sources} ${mcpele_sources} ${PYTHON_LIBS})
add_library(_action_cpp SHARED ${action_source})
target_link_libraries(_action_cpp mcpele_lib)
set_target_properties(_action_cpp PROPERTIES PREFIX "")
