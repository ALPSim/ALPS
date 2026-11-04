foreach(program IN ITEMS convert2xml compactrun snap2vtk)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
install(TARGETS convert2xml snap2vtk
  EXPORT ALPSApplicationTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

install(TARGETS compactrun
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ALPS::convert2xml ALPS::snap2vtk)

if(ALPS_BUILD_TESTING)
  add_test(NAME convert2xml_mc_task
    COMMAND ${CMAKE_COMMAND}
      "-Dconvert2xml=$<TARGET_FILE:convert2xml>"
      "-Dxml_dir=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/mc-task"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/convert2xml_mc_task.cmake")
endif()
