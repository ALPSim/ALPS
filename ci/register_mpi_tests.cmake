# Register retained native MPI fixtures for the focused CI preset. Deferring
# until the root directory finishes makes the normal ALPS targets available.
function(alps_ci_register_mpi_tests)
  if(NOT ALPS_ENABLE_MPI OR NOT ALPS_BUILD_TESTING)
    message(FATAL_ERROR "The MPI CI preset requires MPI and native tests")
  endif()
  foreach(component IN ITEMS alea legacy_parameters)
    if(component STREQUAL "alea")
      set(name observableset_mpi)
    else()
      set(name parameters_mpi)
    endif()
    add_executable(${name} "${CMAKE_SOURCE_DIR}/src/alps/${component}/tests/${name}.C")
    target_link_libraries(${name} PRIVATE ALPS::alps)
  endforeach()
  foreach(case IN ITEMS "observableset_mpi:2:alea" "parameters_mpi:2:legacy_parameters"
      "collect_mpi:2:parapack" "collect_mpi:3:parapack" "comm_mpi:4:parapack"
      "filelock_mpi:2:parapack" "filelock_mpi:3:parapack"
      "halt_mpi:2:parapack" "halt_mpi:3:parapack"
      "info_test_mpi:2:parapack" "info_test_mpi:3:parapack" "process_mpi:8:parapack")
    string(REPLACE ":" ";" fields "${case}")
    list(GET fields 0 name)
    list(GET fields 1 ranks)
    list(GET fields 2 component)
    set(test "ci.mpi.${name}.${ranks}")
    add_test(NAME "${test}" COMMAND "${CMAKE_COMMAND}"
      "-DEXECUTABLE=$<TARGET_FILE:${name}>" "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
      "-DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}" "-DRANKS=${ranks}"
      "-DPREFLAGS=${MPIEXEC_PREFLAGS}" "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
      "-DSOURCE=${CMAKE_SOURCE_DIR}/src/alps/${component}/tests"
      "-DOUTPUT=${CMAKE_BINARY_DIR}/mpi/${test}" "-DNAME=${name}"
      -P "${CMAKE_SOURCE_DIR}/ci/run_mpi_test.cmake")
    set_tests_properties("${test}" PROPERTIES LABELS mpi PROCESSORS "${ranks}" TIMEOUT 300)
  endforeach()
endfunction()
cmake_language(DEFER CALL alps_ci_register_mpi_tests)
