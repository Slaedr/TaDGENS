set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# libraries to be compiled
add_library(base utilities/adatastructures.cpp)

add_library(mesh mesh/amesh2dh.cpp)
target_link_libraries(mesh base)

add_library(fem fem/aelements.cpp fem/aquadrature.cpp)
target_link_libraries(fem mesh)

add_library(spatial spatial/aoutput.cpp spatial/aspatial.cpp)
target_link_libraries(spatial fem)

add_library(spatial_poisson spatial/aspatialpoisson.cpp)
target_link_libraries(spatial_poisson spatial)

add_library(spatial_advection spatial/aspatialadvection.cpp)
target_link_libraries(spatial_advection spatial)

add_library(solvers solvers/atimesteady.cpp)
target_link_libraries(solvers spatial)

#add_library(time time/atimetvdrk.cpp)
#target_link_libraries(time solvers)

#add_library(tadgens_core)
#target_link_libraries(tadgens_core solvers spatial fem mesh base)

add_executable(grid_conv_steady utilities/grid_conv_steady.cpp)
target_link_libraries(grid_conv_steady solvers spatial_advection spatial_poisson)

add_executable(grid_conv_poissonsip poissonSIP.cpp)
target_link_libraries(grid_conv_poissonsip solvers spatial_poisson)

#add_executable(grid_conv_unsteady utilities/grid_conv_steady.cpp)
#target_link_libraries(grid_conv_unsteady tadgens_core)
