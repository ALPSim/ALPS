foreach(program IN ITEMS lattice2xml printgraph)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
install(TARGETS lattice2xml printgraph
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

# pltgraph.py and the historical local fixtures remain source-only.
alps_register_module(NAME tools_lattice
  TARGETS lattice2xml printgraph
  DEPENDS lattice xml)
