# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

alps_add_test(lattice2xml)
set_property(TEST lattice2xml PROPERTY LABELS "cli;compatibility")

foreach(app sparsediag worm)
  if(NOT TARGET ${app})
    continue()
  endif()
  add_test(NAME scheduler_task_failure_${app}
    COMMAND "${CMAKE_COMMAND}" "-DAPP=$<TARGET_FILE:${app}>"
      "-DXML_PATH=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/scheduler_task_failure_${app}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/scheduler_task_failure.cmake")
  set_tests_properties(scheduler_task_failure_${app} PROPERTIES LABELS "cli;scheduler"
    TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

find_package(Python 3.10 QUIET COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  set(archive_arguments "$<TARGET_FILE:txt2archive>" "$<TARGET_FILE:xml2archive>")
  if(TARGET archive)
    list(APPEND archive_arguments "$<TARGET_FILE:archive>")
  endif()
  add_test(NAME archive_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/archive_contract.py"
      ${archive_arguments})
  set_tests_properties(archive_cli_contract PROPERTIES TIMEOUT 60 LABELS "cli;compatibility")
  if(TARGET p2h5)
    add_test(NAME p2h5_cli
      COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_p2h5.py" "$<TARGET_FILE:p2h5>")
    set_tests_properties(p2h5_cli PROPERTIES TIMEOUT 60 LABELS "cli")
  endif()
  if(TARGET spinmc)
    add_test(NAME scheduler_sweeps_spinmc
      COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_scheduler_sweeps.py"
        "$<TARGET_FILE:spinmc>" "${PROJECT_SOURCE_DIR}/src/alps/resources")
    set_tests_properties(scheduler_sweeps_spinmc PROPERTIES TIMEOUT 120 LABELS "cli;scheduler")
    add_test(NAME spinmc_nonsymmetric_coupling
      COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_spinmc_coupling.py"
        "$<TARGET_FILE:spinmc>" "${PROJECT_SOURCE_DIR}/src/alps/resources")
    set_tests_properties(spinmc_nonsymmetric_coupling PROPERTIES TIMEOUT 120 LABELS "cli;compatibility")
  endif()
  if(TARGET qwl AND TARGET qwl_evaluate AND ALPS_ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    add_test(NAME scheduler_qwl_runs_mpi
      COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_scheduler_qwl_runs.py"
        "$<TARGET_FILE:qwl>" "$<TARGET_FILE:qwl_evaluate>" "${PROJECT_SOURCE_DIR}/src/alps/resources"
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "3" ${MPIEXEC_PREFLAGS})
    set_tests_properties(scheduler_qwl_runs_mpi PROPERTIES TIMEOUT 300 PROCESSORS 3 LABELS "cli;scheduler;mpi")
  endif()
endif()
