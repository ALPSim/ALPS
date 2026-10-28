# SPDX-License-Identifier: MIT
#
# Run dmrg twice on a task with the same SEED. The two results must be
# identical; a different SEED changes the last digits of the energy.
# Usage: cmake -DAPP=<dmrg> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

file(REMOVE_RECURSE "${WORK_DIR}")
set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)

foreach(run 1 2)
  file(MAKE_DIRECTORY "${WORK_DIR}/${run}")
  file(WRITE "${WORK_DIR}/${run}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">open chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">8</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"SWEEPS\">2</PARAMETER>
<PARAMETER name=\"MAXSTATES\">20</PARAMETER>
<PARAMETER name=\"SEED\">12345</PARAMETER>
</PARAMETERS></SIMULATION>
")
  file(WRITE "${WORK_DIR}/${run}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
    "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")
  execute_process(COMMAND "${APP}" --write-xml job.in.xml WORKING_DIRECTORY "${WORK_DIR}/${run}"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "run ${run}: exit status ${result}\nstdout:\n${out}\nstderr:\n${err}")
  endif()
  file(READ "${WORK_DIR}/${run}/task.out.xml" result_${run})
endforeach()

if(NOT result_1 STREQUAL result_2)
  message(FATAL_ERROR "two runs with the same SEED differ\nrun 1:\n${result_1}\nrun 2:\n${result_2}")
endif()
