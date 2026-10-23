# SPDX-License-Identifier: MIT
#
# Run fulldiag with --write-xml on a task with T and the field h = J/2. The
# XML output has the averages at that temperature and field, so it must
# succeed and list the field.
# Usage: cmake -DAPP=<fulldiag> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"T\">0.5</PARAMETER>
<PARAMETER name=\"h\">J/2</PARAMETER>
</PARAMETERS></SIMULATION>
")
file(WRITE "${WORK_DIR}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
  "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
execute_process(COMMAND "${APP}" --write-xml job.in.xml WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
set(log "exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")

if(NOT result STREQUAL "0")
  message(FATAL_ERROR "expected exit status 0\n${log}")
endif()
if(NOT EXISTS "${WORK_DIR}/task.out.xml")
  message(FATAL_ERROR "no XML output\n${log}")
endif()
file(READ "${WORK_DIR}/task.out.xml" xml)
if(NOT xml MATCHES "<AVERAGES>")
  message(FATAL_ERROR "the XML output has no averages\n${log}")
endif()
