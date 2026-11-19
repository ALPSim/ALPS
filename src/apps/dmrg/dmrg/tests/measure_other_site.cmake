# SPDX-License-Identifier: MIT
#
# MEASURE_LOCAL and MEASURE_AVERAGE take a site expression in i or the name of
# a bond operator. An expression on another site, such as exchange(i,j), used
# to crash dmrg without a message (issue #385). dmrg must stop with an error
# that names the expression.
# Usage: cmake -DAPP=<dmrg> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)

foreach(measurement "MEASURE_LOCAL[x]|exchange(i,j)" "MEASURE_AVERAGE[x]|Sz(i)*Sz(j)")
  string(REPLACE "|" ";" measurement "${measurement}")
  list(GET measurement 0 name)
  list(GET measurement 1 expression)
  file(REMOVE_RECURSE "${WORK_DIR}")
  file(MAKE_DIRECTORY "${WORK_DIR}")
  file(WRITE "${WORK_DIR}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">open chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"local_S\">1/2</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"CONSERVED_QUANTUMNUMBERS\">Sz</PARAMETER>
<PARAMETER name=\"Sz_total\">0</PARAMETER>
<PARAMETER name=\"SWEEPS\">2</PARAMETER>
<PARAMETER name=\"MAXSTATES\">20</PARAMETER>
<PARAMETER name=\"${name}\">${expression}</PARAMETER>
</PARAMETERS></SIMULATION>
")
  file(WRITE "${WORK_DIR}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
    "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")

  execute_process(COMMAND "${APP}" job.in.xml WORKING_DIRECTORY "${WORK_DIR}"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(log "${name}=${expression}\nexit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")

  if(result STREQUAL "0")
    message(FATAL_ERROR "expected a non-zero exit status\n${log}")
  endif()
  string(FIND "${err}" "'${expression}'" at)
  if(at EQUAL -1 OR NOT err MATCHES "bond operator")
    message(FATAL_ERROR "stderr does not name the expression and bond operators\n${log}")
  endif()
endforeach()
