# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

alps_add_test(lattice2xml)

if(ALPS_CITATION_PYTHON)
add_executable(citation_runtime citation_runtime.cpp)
target_link_libraries(citation_runtime PRIVATE ALPS::alps)
add_test(NAME citation_runtime_serial
  COMMAND "${ALPS_CITATION_PYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/check_citation_runtime.py"
          "$<TARGET_FILE:citation_runtime>")
set_property(TEST citation_runtime_serial PROPERTY LABELS citations)
if(ALPS_HAVE_MPI AND MPIEXEC_EXECUTABLE)
  add_test(NAME citation_runtime_mpi
    COMMAND "${ALPS_CITATION_PYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/check_citation_runtime.py"
            "$<TARGET_FILE:citation_runtime>"
            "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2" ${MPIEXEC_PREFLAGS})
  set_tests_properties(citation_runtime_mpi PROPERTIES LABELS citations PROCESSORS 2)
endif()

endif()

find_package(Python 3.11 QUIET COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  set(archive_arguments "$<TARGET_FILE:txt2archive>" "$<TARGET_FILE:xml2archive>")
  if(TARGET archive)
    list(APPEND archive_arguments "$<TARGET_FILE:archive>")
  endif()
  add_test(NAME archive_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/archive_contract.py"
      ${archive_arguments})
  set_tests_properties(archive_cli_contract PROPERTIES TIMEOUT 60 LABELS "cli;compatibility")
endif()
