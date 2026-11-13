# SPDX-License-Identifier: MIT
#
# Run dmrg on the same task with 1 and with 4 OpenMP threads. The product
# terms are applied in parallel, one thread per output block, in the same
# order as the serial loop, so the two results must be identical. The
# next-nearest-neighbour coupling gives terms between non-adjacent blocks.
# The task is kept small: a BLAS built with OpenMP (such as Homebrew's
# OpenBLAS) follows OMP_NUM_THREADS too, and on large matrices its threaded
# sums would change the last digits on their own.
# Usage: cmake -DAPP=<dmrg> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

file(REMOVE_RECURSE "${WORK_DIR}")
set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
set(ENV{OPENBLAS_NUM_THREADS} 1)

foreach(threads 1 4)
  file(MAKE_DIRECTORY "${WORK_DIR}/${threads}")
  file(WRITE "${WORK_DIR}/${threads}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">nnn open chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">12</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"J'\">0.5</PARAMETER>
<PARAMETER name=\"SWEEPS\">2</PARAMETER>
<PARAMETER name=\"MAXSTATES\">40</PARAMETER>
<PARAMETER name=\"SEED\">12345</PARAMETER>
</PARAMETERS></SIMULATION>
")
  file(WRITE "${WORK_DIR}/${threads}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
    "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")
  set(ENV{OMP_NUM_THREADS} ${threads})
  execute_process(COMMAND "${APP}" --write-xml job.in.xml WORKING_DIRECTORY "${WORK_DIR}/${threads}"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "${threads} threads: exit status ${result}\nstdout:\n${out}\nstderr:\n${err}")
  endif()
  file(READ "${WORK_DIR}/${threads}/task.out.xml" result_${threads})
endforeach()

if(NOT result_1 STREQUAL result_4)
  message(FATAL_ERROR "1 and 4 OpenMP threads give different results\n1 thread:\n${result_1}\n4 threads:\n${result_4}")
endif()
