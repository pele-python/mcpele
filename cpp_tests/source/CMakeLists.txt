cmake_minimum_required(VERSION 2.8)

# Export compile commands for clang tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mcpele_test)

enable_language(CXX)

# TODO: Add a build type option specifically for memory leak testing
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3 -O0")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb3 -O0")
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL Coverage)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -fprofile-arcs -ftest-coverage")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG -fprofile-arcs -ftest-coverage")
endif()

# cmake_policy(SET CMP0015 NEW)

# Add and compile the gtest library
set(gtest_src ../../../pele/cpp_tests/gtest)
set(gtest_binary_dir googletest)
set(gtest_include ${gtest_src}/include/)
include_directories(${gtest_include})
message("gtest include directory ${gtest_include}")
add_subdirectory(../../../pele/cpp_tests/gtest ${gtest_binary_dir})
set(gtest_libs gtest gtest_main)
message("gtest_libs ${gtest_libs}")

# set the mcpele include directory
set(mcpele_include ../../source)
include_directories(${mcpele_include})
message("mcpele include directory ${mcpele_include}")

# build the mcpele library
FILE(GLOB mcpele_sources ${mcpele_include}/*.cpp)
add_library(mcpele_lib ${mcpele_sources})

set(PELE_DIR ../../../pele)

# set the pele include directory
set(pele_include ${PELE_DIR}/source CACHE STRING "the pele c++ source directory")

if(EXISTS ${pele_include}/pele/array.hpp)
	message("pele include directory: ${pele_include}")
else()
	message(FATAL_ERROR "pele include directory is not correct: ${pele_include} : use ccmake to set it")
endif()

set(PELE_EXTERN_DIR ${PELE_DIR}/extern/install)
set(EXTERN_INCLUDE_DIR ${PELE_EXTERN_DIR}/include/) # FOR EIGEN, CVODE etc
set(LAPACKE_INCLUDE_DIR ${PELE_EXTERN_DIR}/include/Eigen/src/misc)

include_directories(${pele_include})
include_directories(${EXTERN_INCLUDE_DIR})
include_directories(${LAPACKE_INCLUDE_DIR})

# build the pele library
FILE(GLOB pele_sources ${pele_include}/*.c*)
add_library(pele_lib ${pele_sources})

# get all the source files
FILE(GLOB sources *.c*)

add_executable(test_main ${sources})
target_link_libraries(test_main pele_lib mcpele_lib gtest gtest_main pthread)

# Custom version of this function which came from CodeCoverage.cmake
# It removes coverage information for the pele library.
FUNCTION(SETUP_TARGET_FOR_COVERAGE_CUSTOM _targetname _testrunner _outputname)
	IF(NOT LCOV_PATH)
		MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
	ENDIF() # NOT LCOV_PATH

	IF(NOT GENHTML_PATH)
		MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
	ENDIF() # NOT GENHTML_PATH

	# Setup target
	ADD_CUSTOM_TARGET(${_targetname}

		# Cleanup lcov
		${LCOV_PATH} --directory . --zerocounters

		# Run tests
		COMMAND ${_testrunner} ${ARGV3}

		# Capturing lcov counters and generating report
		COMMAND ${LCOV_PATH} --directory . --capture --output-file ${_outputname}.info
		COMMAND ${LCOV_PATH} --remove ${_outputname}.info 'tests/*' '/usr/*' '/pele/*' --output-file ${_outputname}.info.cleaned
		COMMAND ${GENHTML_PATH} -o ${_outputname} ${_outputname}.info.cleaned
		COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned

		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
	)

	# Show info where to find the report
	ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Open ./${_outputname}/index.html in your browser to view the coverage report."
	)
ENDFUNCTION() # SETUP_TARGET_FOR_COVERAGE_CUSTOM

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
	include(CodeCoverage)
	setup_target_for_coverage_custom(mcpele_coverage test_main coverage)
endif()
