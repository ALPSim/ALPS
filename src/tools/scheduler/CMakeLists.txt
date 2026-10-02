foreach(program IN ITEMS convert2xml compactrun snap2vtk)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
install(TARGETS convert2xml compactrun snap2vtk
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

alps_register_module(NAME tools_scheduler
  TARGETS convert2xml compactrun snap2vtk
  DEPENDS osiris scheduler utilities xml)
