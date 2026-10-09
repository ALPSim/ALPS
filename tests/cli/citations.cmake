# Citation integration checks are optional for native-only builds.
if(NOT ALPS_CITATION_PYTHON)
  find_package(Python3 3.11 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    set(ALPS_CITATION_PYTHON "${Python3_EXECUTABLE}")
  endif()
endif()
if(ALPS_CITATION_PYTHON)
  execute_process(COMMAND "${ALPS_CITATION_PYTHON}" -c "import yaml, jsonschema"
    RESULT_VARIABLE missing_dependencies OUTPUT_QUIET ERROR_QUIET)
  if(NOT missing_dependencies EQUAL 0)
    set(ALPS_CITATION_PYTHON "")
  endif()
endif()
if(NOT ALPS_CITATION_PYTHON)
  message(STATUS "Skipping citation integration checks: Python with PyYAML/jsonschema unavailable")
  return()
endif()

add_test(NAME citation_catalog COMMAND "${ALPS_CITATION_PYTHON}"
  "${PROJECT_SOURCE_DIR}/.github/scripts/generate_citations.py" --check)
set_tests_properties(citation_catalog PROPERTIES LABELS citations TIMEOUT 60)

add_executable(citation_runtime "${CMAKE_CURRENT_LIST_DIR}/citation_runtime.cpp")
target_link_libraries(citation_runtime PRIVATE ALPS::alps)
add_test(NAME citation_runtime_serial COMMAND "${ALPS_CITATION_PYTHON}"
  "${CMAKE_CURRENT_LIST_DIR}/check_citation_runtime.py" "$<TARGET_FILE:citation_runtime>")
set_tests_properties(citation_runtime_serial PROPERTIES LABELS citations TIMEOUT 180)
if(ALPS_ENABLE_MPI AND MPIEXEC_EXECUTABLE)
  add_test(NAME citation_runtime_mpi COMMAND "${ALPS_CITATION_PYTHON}"
    "${CMAKE_CURRENT_LIST_DIR}/check_citation_runtime.py" "$<TARGET_FILE:citation_runtime>"
    "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2" ${MPIEXEC_PREFLAGS})
  set_tests_properties(citation_runtime_mpi PROPERTIES LABELS "citations;mpi" PROCESSORS 2 TIMEOUT 180)
endif()
foreach(app dmrg qwl dirloop_sse worm spinmc fulldiag sparsediag hirschfye loop
    dmft hybridization interaction fulldiag_evaluate spinmc_evaluate checksign
    qwl_evaluate worm_evaluate simplemc pevaluate maxent)
  if(NOT TARGET ${app})
    continue()
  endif()
  string(REPLACE "_evaluate" "" component "${app}")
  if(app STREQUAL "loop")
    set(component looper)
  elseif(app MATCHES "^(checksign|simplemc|pevaluate|maxent)$")
    set(component framework)
  endif()
  add_test(NAME citation_cli_${app} COMMAND "${ALPS_CITATION_PYTHON}"
    "${CMAKE_CURRENT_LIST_DIR}/check_citation_notice.py" "$<TARGET_FILE:${app}>" "${component}")
  set_tests_properties(citation_cli_${app} PROPERTIES LABELS "citations;cli" TIMEOUT 300)
  if(ALPS_ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    add_test(NAME citation_cli_mpi_${app} COMMAND "${ALPS_CITATION_PYTHON}"
      "${CMAKE_CURRENT_LIST_DIR}/check_citation_notice.py" "$<TARGET_FILE:${app}>" "${component}"
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2" ${MPIEXEC_PREFLAGS})
    set_tests_properties(citation_cli_mpi_${app} PROPERTIES LABELS "citations;cli;mpi" PROCESSORS 2 TIMEOUT 300)
  endif()
endforeach()
if(TARGET loop)
  add_test(NAME citation_calculation_loop COMMAND "${ALPS_CITATION_PYTHON}"
    "${CMAKE_CURRENT_LIST_DIR}/check_citation_startup.py" "$<TARGET_FILE:loop>"
    "${PROJECT_SOURCE_DIR}/src/alps/resources")
  set_tests_properties(citation_calculation_loop PROPERTIES LABELS "citations;cli" TIMEOUT 120)
endif()
if(TARGET dmft)
  add_test(NAME citation_dmft_builtin_solvers COMMAND "${ALPS_CITATION_PYTHON}"
    "${CMAKE_CURRENT_LIST_DIR}/check_citation_dmft.py" "$<TARGET_FILE:dmft>")
  set_tests_properties(citation_dmft_builtin_solvers PROPERTIES LABELS "citations;cli" TIMEOUT 120)
endif()
