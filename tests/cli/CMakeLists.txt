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
