add_executable(testelements testelementtri.cpp)
target_link_libraries(testelements fem mesh base)

configure_file(../common_inputs/trimesh-skew_p2.msh trimesh-skew_p2.msh)

# add_test(NAME P2_QuadraturePoints
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testelements
#   )
