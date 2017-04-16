cmake_minimum_required(VERSION 2.8)

project (tadgens)

# to compile with OpenMP
if(OMP)
	message(STATUS "Compiling with OpenMP")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
	else()
		message(SEND_ERROR "Don't know OpenMP flag for this compiler")
	endif()
endif()

# set compile options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftree-vectorizer-verbose=2")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Eigen
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE")
include_directories($ENV{EIGEN_DIR})

# ---------------------------------------------------------------------------- #

# flags and stuff

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ggdb -DDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -UDEBUG -DNDEBUG")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Debug mode")
else()
	message(STATUS "Release mode")
endif()

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# libraries to be compiled
add_library(tadgens_base aoutput.cpp aspatial.cpp aelements.cpp aquadrature.cpp amesh2dh.cpp adatastructures.cpp)

#add_library(tadgens_poisson aspatialpoisson.cpp)
#target_link_libraries(tadgens_poisson tadgens_base)

add_library(tadgens_advection aspatialadvection.cpp)
target_link_libraries(tadgens_advection tadgens_base)

# for the final executable(s)

#add_subdirectory(utilities)

#add_executable(poissonsip poissonSIP.cpp)
#target_link_libraries(poissonsip tadgens_poisson)
add_executable(advect ladvection.cpp)
target_link_libraries(advect tadgens_advection)
#add_executable(poissonc poissonC.cpp)
#target_link_libraries(poissonc tadgens_poisson)

