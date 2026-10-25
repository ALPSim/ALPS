# SPDX-License-Identifier: MIT
#
# Run fulldiag with --write-xml, and then fulldiag_evaluate, on a task with T
# and the field h = J/2, and on the same task with h = 0.5. Both must succeed,
# and the averages and the evaluated plots must agree. The plots are over a
# range of fields measured from h, so a wrong h would shift them. A blank or
# unevaluable temperature or field must stop fulldiag_evaluate.
# Usage: cmake -DAPP=<fulldiag> -DEVALUATE=<fulldiag_evaluate> -DXML_PATH=<resources>
#              -DWORK_DIR=<dir> -P ...

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
file(REMOVE_RECURSE "${WORK_DIR}")

# runs the task with field h in WORK_DIR/<name>; sets <name>_averages to the
# averages of the XML output and <name>_plot to the evaluated magnetization
function(run_task name h)
  set(dir "${WORK_DIR}/${name}")
  file(MAKE_DIRECTORY "${dir}")
  file(WRITE "${dir}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"CONSERVED_QUANTUMNUMBERS\">Sz</PARAMETER>
<PARAMETER name=\"T\">0.5</PARAMETER>
<PARAMETER name=\"h\">${h}</PARAMETER>
</PARAMETERS></SIMULATION>
")
  file(WRITE "${dir}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
    "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")

  execute_process(COMMAND "${APP}" --write-xml job.in.xml WORKING_DIRECTORY "${dir}"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(log "h = ${h}: fulldiag exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "expected exit status 0\n${log}")
  endif()
  if(NOT EXISTS "${dir}/task.out.xml")
    message(FATAL_ERROR "no XML output\n${log}")
  endif()
  file(READ "${dir}/task.out.xml" xml)
  if(NOT xml MATCHES "<AVERAGES>.*</AVERAGES>")
    message(FATAL_ERROR "the XML output has no averages\n${log}")
  endif()
  set(${name}_averages "${CMAKE_MATCH_0}" PARENT_SCOPE)

  execute_process(COMMAND "${EVALUATE}" --T_MIN 0.1 --T_MAX 1 --DELTA_T 0.1
    --H_MIN 0 --H_MAX 1 --DELTA_H 0.5 task.out.xml
    WORKING_DIRECTORY "${dir}" RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(log "h = ${h}: fulldiag_evaluate exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "expected exit status 0\n${log}")
  endif()
  if(NOT EXISTS "${dir}/task.plot.magnetization.xml")
    message(FATAL_ERROR "no magnetization plot\n${log}")
  endif()
  file(READ "${dir}/task.plot.magnetization.xml" plot)
  set(${name}_plot "${plot}" PARENT_SCOPE)
endfunction()

run_task(expression "J/2")
run_task(number "0.5")

if(NOT expression_averages STREQUAL number_averages)
  message(FATAL_ERROR "the averages differ\nh = J/2:\n${expression_averages}\nh = 0.5:\n${number_averages}")
endif()
if(NOT expression_averages MATCHES "Magnetization")
  message(FATAL_ERROR "the averages have no magnetization\n${expression_averages}")
endif()
if(NOT expression_plot STREQUAL number_plot)
  message(FATAL_ERROR "the evaluated magnetizations differ\nh = J/2:\n${expression_plot}\nh = 0.5:\n${number_plot}")
endif()

# a temperature or field that is blank or cannot be evaluated is an error
# naming it, not a plot of NaN
foreach(case "T; " "h;K/2")
  list(GET case 0 name)
  list(GET case 1 value)
  execute_process(COMMAND "${EVALUATE}" --T 0.5 --${name} "${value}" task.out.xml
    WORKING_DIRECTORY "${WORK_DIR}/expression" RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(result STREQUAL "0" OR NOT err MATCHES "parameter '${name}'")
    message(FATAL_ERROR "${name} = '${value}': expected an error naming the parameter\n"
      "exit status: ${result}\nstderr:\n${err}")
  endif()
endforeach()
