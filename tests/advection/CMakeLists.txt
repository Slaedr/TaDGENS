configure_file(advect-t.control advect-t.control)
configure_file(advect-l.control advect-l.control)
configure_file(advect-l-quad.control advect-l-quad.control)
configure_file(advect-t-struct.control advect-t-struct.control)

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
  
  set(steadyadvection_quad_meshfiles
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad0.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad1.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad2.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad3.msh 
	#${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad4.msh 
	)
  
set(steadyadvection_struct_meshfiles
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct0.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct1.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct2.msh 
	${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct3.msh 
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
  
  foreach(imesh 0 1 2 3 4)
	add_custom_command(
	  COMMENT "Build quad grids needed for steady advection test"
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad${imesh}.msh
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2
	    -o ${CMAKE_CURRENT_BINARY_DIR}/grids/squarequad${imesh}.msh
	    -v 2
	    square-quad.geo
	  MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/square-quad.geo"
	  )
  endforeach(imesh)
  
  foreach(imesh 0 1 2 3 4)
	add_custom_command(
	  COMMENT "Build structured grids needed for steady advection test"
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct${imesh}.msh
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  COMMAND ${GMSH_EXEC} -setnumber ref ${imesh} -2 -format msh2
	    -o ${CMAKE_CURRENT_BINARY_DIR}/grids/squarestruct${imesh}.msh
	    -v 2
	    square-struct.geo
	  MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/square-struct.geo"
	  )
  endforeach(imesh)

  add_custom_target(steadyadvection_meshes DEPENDS ${steadyadvection_tri_meshfiles} 
                                                   ${steadyadvection_quad_meshfiles}
												   ${steadyadvection_struct_meshfiles})
  add_dependencies(grid_conv_steady steadyadvection_meshes)

  add_test(NAME SteadyAdvection_SolutionConvergence_Lagrange_P1
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
	${CMAKE_CURRENT_BINARY_DIR}/advect-l.control
	)
  
  add_test(NAME SteadyAdvection_SolutionConvergence_Lagrange_P1_Quad
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
	${CMAKE_CURRENT_BINARY_DIR}/advect-l-quad.control
	)
  
  add_test(NAME SteadyAdvection_SolutionConvergence_Taylor_P1
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
  	${CMAKE_CURRENT_BINARY_DIR}/advect-t.control
	)
  
add_test(NAME SteadyAdvection_SolutionConvergence_Taylor_P1_Struct
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  	COMMAND ${SEQEXEC} ${SEQTASKS} ${CMAKE_BINARY_DIR}/grid_conv_steady
  	${CMAKE_CURRENT_BINARY_DIR}/advect-t-struct.control
	)

endif()
