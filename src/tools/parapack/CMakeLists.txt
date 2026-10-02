foreach(program IN ITEMS pevaluate poutput)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
install(TARGETS pevaluate poutput
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

# xml2archive remains unbuilt; inventorying it does not add a target.
alps_register_module(NAME tools_parapack
  TARGETS pevaluate poutput
  EXTRA_FILES xml2archive.C
  DEPENDS parapack)
