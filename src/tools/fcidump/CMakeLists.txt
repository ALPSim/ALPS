# Export of ALPS fermion models to the FCIDUMP integral format.
add_executable(alps2fcidump alps2fcidump.C)
target_link_libraries(alps2fcidump PRIVATE ALPS::alps)
install(TARGETS alps2fcidump
  EXPORT ALPSApplicationTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ALPS::alps2fcidump)

if(ALPS_BUILD_TESTING)
  add_test(NAME alps2fcidump_help COMMAND alps2fcidump --help)
  set_tests_properties(alps2fcidump_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:"
    LABELS "cli" TIMEOUT 60)
endif()
