
add_executable(poissonc poissonC.cpp aspatialpoisson_continuous.cpp)
target_link_libraries(poissonc spatial)

configure_file(poisson-continuous-tri.control poisson-continuous-tri.control)

set(poisson_tri_meshfiles
  ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri0.msh 
  ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri1.msh 
  ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri2.msh 
  ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri3.msh 
  ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri4.msh 
  )

foreach(imesh 0 1 2 3 4)
  add_custom_command(
	COMMENT "Build grid needed for Poisson test"
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri${imesh}.msh
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/grids
	COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2
	-o ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri${imesh}.msh
	-v 2
	grids/square-tri.geo
	MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/grids/square-tri.geo"
	)
endforeach(imesh)

add_custom_target(poisson_meshes DEPENDS ${poisson_tri_meshfiles})

add_dependencies(poissonc poisson_meshes)

add_test(NAME Poisson_C_SolutionConvergence_Lagrange_P1_Tri
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_CURRENT_BINARY_DIR}/poissonc
  ${CMAKE_CURRENT_BINARY_DIR}/poisson-continuous-tri.control
  )
