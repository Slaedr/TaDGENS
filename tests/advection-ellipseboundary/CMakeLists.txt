configure_file(advect-l-p1.control advect-l-p1.control)
configure_file(advect-l-p2.control advect-l-p2.control)

if(${GMSH_EXEC} STREQUAL "GMSH_EXEC-NOTFOUND")
  message(WARNING "Steady advection test not built because Gmsh was not found")
else()

  set(steadyadvectioncurved_tri_meshfiles
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_0.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_1.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_2.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_3.msh 
	#${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_4.msh 
	)
  set(steadyadvectioncurved_p2_tri_meshfiles
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_0.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_1.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_2.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_3.msh 
	#${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_4.msh 
	)

  foreach(imesh 0 1 2 3 4)
	add_custom_command(
	  COMMENT "Build p1 grids needed for steady advection curved outflow test"
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_${imesh}.msh
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/grids
	  COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2
	    -o ${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_${imesh}.msh
	    -v 2
	    square-ellipse.geo
	  MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/square-ellipse.geo"
	  )
	add_custom_command(
	  COMMENT "Build p2 grids needed for steady advection curved outflow test"
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_${imesh}.msh
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/grids
	  COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2 -order 2
	    -o ${CMAKE_CURRENT_BINARY_DIR}/grids/squareellipse_p2_${imesh}.msh
	    -v 2
	    square-ellipse.geo
	  MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/square-ellipse.geo"
	  )
  endforeach(imesh)

  add_custom_target(steadyadvectioncurved_meshes DEPENDS ${steadyadvectioncurved_tri_meshfiles})
  add_dependencies(grid_conv_steady steadyadvectioncurved_meshes)
  add_custom_target(steadyadvectioncurved_p2_meshes DEPENDS ${steadyadvectioncurved_p2_tri_meshfiles})
  add_dependencies(grid_conv_steady steadyadvectioncurved_p2_meshes)

  add_test(NAME SteadyAdvectionEllipseOutflow_SolutionConvergence_Lagrange_P1
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
	${CMAKE_CURRENT_BINARY_DIR}/advect-l-p1.control
	)
  add_test(NAME SteadyAdvectionEllipseOutflow_SolutionConvergence_Lagrange_P2
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
	${CMAKE_CURRENT_BINARY_DIR}/advect-l-p2.control
	)

endif()
