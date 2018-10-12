cmake_minimum_required(VERSION 2.8)

project (tadgens)

# to compile with OpenMP
if(NOOMP)
  message(STATUS "Compiling without OpenMP")
else()
  message(STATUS "Compiling with OpenMP")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  else()
	message(SEND_ERROR "Don't know OpenMP flag for this compiler")
  endif()
endif()

# set common compile options
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Eigen
find_path(EIGEN3_INCLUDE_DIR signature_of_eigen3_matrix_library
  PATHS ENV EIGEN3_ROOT DOC "Eigen3 include directory")
if(${EIGEN3_INCLUDE_DIR} STREQUAL "EIGEN3_INCLUDE_DIR-NOTFOUND")
  message(SEND_ERROR "Eigen3 not found!")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "Found Eigen 3 at ${EIGEN3_INCLUDE_DIR}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE")


# ---------------------------------------------------------------------------- #

# flags and stuff

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -DDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -UDEBUG -DNDEBUG")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Debug mode")
else()
	message(STATUS "Release mode")
endif()

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# for the final executable(s)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
