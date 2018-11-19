find_program(GMSH_EXEC gmsh)
if(${GMSH_EXEC} STREQUAL "GMSH_EXEC-NOTFOUND")
  message(WARNING "Gmsh not found")
else()
  message(STATUS "Found Gmsh at ${GMSH_EXEC}")
endif()

set(SEQEXEC "")
set(SEQTASKS "")

add_subdirectory(fem)
add_subdirectory(advection)
add_subdirectory(advection-ellipseboundary)
add_subdirectory(poisson)
