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

  set(fcidump_tests "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  set(fcidump_work "${CMAKE_CURRENT_BINARY_DIR}/test-work")
  file(MAKE_DIRECTORY "${fcidump_work}")
  set(fcidump_xml "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")

  # A periodic Hubbard chain with hopping, U, V and mu against a checked-in FCIDUMP.
  add_test(NAME alps2fcidump_hubbard
    COMMAND "${CMAKE_COMMAND}" -DTOOL=$<TARGET_FILE:alps2fcidump>
      -DINPUT=${fcidump_tests}/hubbard.par -DREFERENCE=${fcidump_tests}/hubbard.fcidump
      -P "${fcidump_tests}/compare.cmake"
    WORKING_DIRECTORY "${fcidump_work}")
  set_tests_properties(alps2fcidump_hubbard PROPERTIES
    ENVIRONMENT_MODIFICATION "${fcidump_xml}" LABELS "cli" TIMEOUT 60)

  # Sector parameters: N_total and Sz_total of "alternative fermion" give NELEC
  # and MS2, and TOTAL_MOMENTUM is ignored.
  foreach(case IN ITEMS
      "constrained_n|NELEC=4, MS2=0"
      "constrained_sz|NELEC=4, MS2=0"
      "momentum|NELEC=4, MS2=0")
    string(REPLACE "|" ";" case "${case}")
    list(GET case 0 name)
    list(GET case 1 expected)
    add_test(NAME alps2fcidump_${name}
      COMMAND alps2fcidump ${name}.par "${fcidump_work}/${name}.fcidump"
      WORKING_DIRECTORY "${fcidump_tests}")
    set_tests_properties(alps2fcidump_${name} PROPERTIES
      ENVIRONMENT_MODIFICATION "${fcidump_xml}" PASS_REGULAR_EXPRESSION "Wrote [^\n]*${expected}"
      LABELS "cli" TIMEOUT 60)
  endforeach()

  # Models an FCIDUMP cannot hold, and electron counts that do not fit, are
  # refused with a message, without a file. They run in the source directory
  # for the test model library; the refused FCIDUMP would go to the build tree.
  foreach(case IN ITEMS
      "fixed_n|BASIS also constrains N "
      "unconstrained|does not fix Nup and Ndown"
      "lone_nup|Nup_total is set without Ndown_total"
      "conflicting_n|Nup_total = 2 and Ndown_total = 2 disagree with N_total = 3"
      "overfilled|6 up and 0 down electrons do not fit in 4 orbitals")
    string(REPLACE "|" ";" case "${case}")
    list(GET case 0 name)
    list(GET case 1 expected)
    add_test(NAME alps2fcidump_refuses_${name}
      COMMAND alps2fcidump ${name}.par "${fcidump_work}/refused_${name}.fcidump"
      WORKING_DIRECTORY "${fcidump_tests}")
    set_tests_properties(alps2fcidump_refuses_${name} PROPERTIES
      ENVIRONMENT_MODIFICATION "${fcidump_xml}" PASS_REGULAR_EXPRESSION "${expected}"
      LABELS "cli" TIMEOUT 60)
  endforeach()
endif()
