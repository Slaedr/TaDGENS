
#add_executable(poissonc poissonC.cpp)
#target_link_libraries(poissonc tadgens_poisson)

add_test(NAME Poisson_SIP_SolutionConvergence_Lagrange_P1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_poissonsip
  ${CMAKE_CURRENT_SOURCE_DIR}/poisson.control
  )

