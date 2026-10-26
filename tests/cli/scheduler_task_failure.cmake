# SPDX-License-Identifier: MIT
#
# Run a scheduler application on a job whose first task names an unknown
# lattice and whose second task file is missing. The third task must still
# run, the failures must be reported on stderr, and the exit status must be 2.
# Diagonalization fails while the scheduler reads the task, Monte Carlo when
# it creates the workers, so run one application of each kind.
# Usage: cmake -DAPP=<application> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

function(write_task name lattice)
  file(WRITE "${WORK_DIR}/${name}.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">${lattice}</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"CONSERVED_QUANTUMNUMBERS\">Sz</PARAMETER>
<PARAMETER name=\"Sz_total\">0</PARAMETER>
<PARAMETER name=\"NUMBER_EIGENVALUES\">1</PARAMETER>
<PARAMETER name=\"T\">1</PARAMETER>
<PARAMETER name=\"THERMALIZATION\">20</PARAMETER>
<PARAMETER name=\"SWEEPS\">200</PARAMETER>
</PARAMETERS></SIMULATION>
")
endfunction()

write_task(bad "no such lattice")
write_task(good "chain lattice")
set(job "<JOB><OUTPUT file=\"job.out.xml\"/>")
foreach(task bad missing good)
  string(APPEND job "<TASK status=\"new\"><INPUT file=\"${task}.in.xml\"/>"
                    "<OUTPUT file=\"${task}.out.xml\"/></TASK>")
endforeach()
file(WRITE "${WORK_DIR}/job.in.xml" "${job}</JOB>\n")

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
execute_process(COMMAND "${APP}" job.in.xml WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
set(log "exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")

if(NOT result STREQUAL "2")
  message(FATAL_ERROR "expected exit status 2\n${log}")
endif()
foreach(expected
    "cannot create task 1 from [^\n]*bad.in.xml: [^\n]*no such lattice"
    "cannot create task 2 from [^\n]*missing.in.xml: file does not exist"
    "2 of 3 tasks could not be created")
  if(NOT err MATCHES "${expected}")
    message(FATAL_ERROR "stderr does not match '${expected}'\n${log}")
  endif()
endforeach()
if(NOT EXISTS "${WORK_DIR}/good.out.h5")
  message(FATAL_ERROR "the valid task did not run\n${log}")
endif()
# A restart from the checkpointed job retries the failed tasks.
file(READ "${WORK_DIR}/job.out.xml" checkpoint)
string(REGEX MATCHALL "status=\"[a-z]+\"" statuses "${checkpoint}")
if(NOT statuses STREQUAL "status=\"new\";status=\"new\";status=\"finished\"")
  message(FATAL_ERROR "unexpected task statuses ${statuses} in\n${checkpoint}")
endif()
