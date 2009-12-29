######################################################################
# applications
######################################################################

add_definitions(${LAPACK_DEFINITIONS})
add_subdirectory(mc)
add_subdirectory(qmc)
add_subdirectory(diag)
add_subdirectory(dmrg)
  

######################################################################
# GUI's
######################################################################

# add_subdirectory(gui/simulator)

