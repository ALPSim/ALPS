######################################################################
# applications
######################################################################

# This is TRUE if libguide is being linked
IF (REQUIRE_PTHREAD)
  SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
ENDIF (REQUIRE_PTHREAD)

add_definitions(${LAPACK_DEFINITIONS})
add_subdirectory(mc)
add_subdirectory(qmc)
add_subdirectory(diag)
add_subdirectory(dmrg)
  

######################################################################
# GUI's
######################################################################

# add_subdirectory(gui/simulator)

