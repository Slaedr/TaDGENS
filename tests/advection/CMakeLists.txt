configure_file(advect-t.control advect-t.control)
configure_file(advect-l.control advect-l.control)

if(${GMSH_EXEC} STREQUAL "GMSH_EXEC-NOTFOUND")
  message(WARNING "Steady advection test not built because Gmsh was not found")
else()

  set(steadyadvection_tri_meshfiles
	${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri0.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri1.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri2.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri3.msh 
	#${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri4.msh 
	)

  foreach(imesh 0 1 2 3 4)
	add_custom_command(
	  COMMENT "Build grid needed for steady advection test"
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri${imesh}.msh
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/grids
	  COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2
	    -o ${CMAKE_CURRENT_BINARY_DIR}/grids/squaretri${imesh}.msh
	    -v 2
	    square-tri.geo
	  MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/square-tri.geo"
	  )
  endforeach(imesh)

  add_custom_target(steadyadvection_meshes DEPENDS ${steadyadvection_tri_meshfiles})
  add_dependencies(grid_conv_steady steadyadvection_meshes)

  add_test(NAME SteadyAdvection_SolutionConvergence_Lagrange_P2
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
	${CMAKE_CURRENT_BINARY_DIR}/advect-l.control
	)
  
  add_test(NAME SteadyAdvection_SolutionConvergence_Taylor_P2
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
  	${CMAKE_CURRENT_BINARY_DIR}/advect-t.control
	)

endif()
