# SPDX-License-Identifier: MIT
#
# Run a diagonalization application on a task whose NUMBER_EIGENVALUES is not
# an integer. The error must name the parameter on stderr, and the
# application must exit with a non-zero status.
# Usage: cmake -DAPP=<application> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"NUMBER_EIGENVALUES\">2.5</PARAMETER>
</PARAMETERS></SIMULATION>
")
file(WRITE "${WORK_DIR}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
  "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
execute_process(COMMAND "${APP}" job.in.xml WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
set(log "exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")

if(result STREQUAL "0")
  message(FATAL_ERROR "expected a non-zero exit status\n${log}")
endif()
if(NOT err MATCHES "parameter 'NUMBER_EIGENVALUES'")
  message(FATAL_ERROR "stderr does not name NUMBER_EIGENVALUES\n${log}")
endif()
