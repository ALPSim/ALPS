# SPDX-License-Identifier: MIT
set(test_command "${mpiexec}" ${mpiexec_numproc_flag} ${procs}
  ${mpiexec_preflags} "${cmd_path}" ${mpiexec_postflags} ${opt})
include("${CMAKE_CURRENT_LIST_DIR}/run_test.cmake")
