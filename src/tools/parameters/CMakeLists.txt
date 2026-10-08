# Existing parameter grammar and simulation job preparation tools.
foreach(program IN ITEMS parameter2xml parameter2hdf5)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
add_executable(p2h5 p2h5.cpp)
target_link_libraries(p2h5 PRIVATE ALPS::alps)
install(TARGETS parameter2xml
  EXPORT ALPSApplicationTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

install(TARGETS parameter2hdf5 p2h5
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ALPS::parameter2xml)
