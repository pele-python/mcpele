cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
SET(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -std=c++0x")

#cmake_policy(SET CMP0015 NEW)

# set the mcpele include directory
set(mcpele_include ${CMAKE_SOURCE_DIR})
include_directories(${mcpele_include})
message("mcpele include directory ${mcpele_include}")

# build the mcpele library
FILE(GLOB mcpele_sources ${mcpele_include}/mcpele/*.cpp)
add_library(mcpele_lib SHARED ${mcpele_sources})

# set the pele include directory
set(pele_include ${CMAKE_SOURCE_DIR}/../../pele/source CACHE STRING "the pele c++ source directory")
if (EXISTS ${pele_include}/pele/array.h)
  message("pele include directory: ${pele_include}")
else()
  message(FATAL_ERROR "pele include directory is not correct: ${pele_include} : use ccmake to set it")
endif()
include_directories(${pele_include})

# build the pele library
FILE(GLOB pele_sources ${pele_include}/*.cpp)
add_library(pele_lib SHARED ${pele_sources})

set(install_dir ${CMAKE_BINARY_DIR}/dest)
install(TARGETS pele_lib mcpele_lib LIBRARY  DESTINATION ${install_dir})
