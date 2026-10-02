# Existing parameter grammar and simulation job preparation tools.
foreach(program IN ITEMS parameter2xml parameter2hdf5)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
add_executable(p2h5 p2h5.cpp)
target_link_libraries(p2h5 PRIVATE ALPS::alps)
install(TARGETS parameter2xml parameter2hdf5 p2h5
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

alps_register_module(NAME tools_parameters
  TARGETS parameter2xml parameter2hdf5 p2h5
  DEPENDS hdf5 legacy_parameters xml)
