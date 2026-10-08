foreach(program IN ITEMS pevaluate poutput)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
install(TARGETS pevaluate poutput
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
