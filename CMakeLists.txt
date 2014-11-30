# this file is mean to be parsed by python before being a proper CMakeLists.txt file
# in particular, the following strings will be set by python
set(numpy_include /usr/lib/python2.7/dist-packages/numpy/core/include)
set(PYTHON_INCLUDE_DIRS /usr/include/python2.7 /usr/include/python2.7)
set(PYTHON_LDFLAGS "-lpthread -ldl -lutil -lm -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions")
set(pele_include /home/stefano/Work/pele/source)

message("python include dirs: ${PYTHON_INCLUDE_DIRS}")
message("numpy include dirs: ${numpy_include}")
message("pele include dirs: ${pele_include}")
include_directories(${numpy_include})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${pele_include})

cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
SET(CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wextra -pedantic -O3")

#cmake_policy(SET CMP0015 NEW)

# set the mcpele include directory
set(mcpele_include ${CMAKE_SOURCE_DIR}/source)
include_directories(${mcpele_include})
message("mcpele include directory: ${mcpele_include}")

# build the mcpele library
FILE(GLOB pele_sources ${pele_include}/*.cpp)
FILE(GLOB mcpele_sources ${mcpele_include}/*.cpp)
add_library(mcpele_lib SHARED ${mcpele_sources} ${pele_sources})

function(make_cython_lib cython_cxx_source)
  get_filename_component(library_name ${cython_cxx_source} NAME)
  string(REGEX REPLACE ".cxx$" "" library_name ${library_name})
  add_library(${library_name} SHARED ${cython_cxx_source})
  target_link_libraries(${library_name} mcpele_lib)
  set_target_properties(${library_name} PROPERTIES PREFIX "")
  message("making library ${library_name} from source ${cython_cxx_source}")
endfunction(make_cython_lib)



make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_pele_mc.cxx)
make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_monte_carlo_cpp.cxx)
make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_takestep_cpp.cxx)
make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_accept_test_cpp.cxx)
make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_conf_test_cpp.cxx)
make_cython_lib(${CMAKE_SOURCE_DIR}/mcpele/monte_carlo/_action_cpp.cxx)
