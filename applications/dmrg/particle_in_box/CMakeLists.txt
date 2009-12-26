add_executable(simple_dmrg dmrgmain.C dmrg3.C)
target_link_libraries(simple_dmrg alps boost ${LAPACK_LIBRARIES})
install(TARGETS simple_dmrg RUNTIME DESTINATION bin)
