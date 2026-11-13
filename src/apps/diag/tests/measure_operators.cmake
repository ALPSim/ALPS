# SPDX-License-Identifier: MIT
#
# Run fulldiag --write-xml with one MEASURE_AVERAGE, MEASURE_LOCAL and
# MEASURE_STRUCTURE_FACTOR key and two MEASURE_CORRELATIONS keys, on the
# Sz_total = 1 sector of the open four-site Heisenberg chain. Every measurement
# must appear under its bracketed name with the exact ground-state value.
# "Splus:Sminus" is not symmetric in its operators: its on-site value is
# 1/2 + <Sz>, and would be 1/2 - <Sz> with the two swapped. The chain is open
# because local measurements are skipped on translation-invariant lattices.
# Usage: cmake -DAPP=<fulldiag> -DXML_PATH=<resources> -DWORK_DIR=<dir> -P ...

set(ENV{ALPS_XML_PATH} "${XML_PATH}")
set(ENV{ALPS_NO_CITATIONS} 1)
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

file(WRITE "${WORK_DIR}/task.in.xml" "<SIMULATION><PARAMETERS>
<PARAMETER name=\"LATTICE\">open chain lattice</PARAMETER>
<PARAMETER name=\"MODEL\">spin</PARAMETER>
<PARAMETER name=\"L\">4</PARAMETER>
<PARAMETER name=\"J\">1</PARAMETER>
<PARAMETER name=\"CONSERVED_QUANTUMNUMBERS\">Sz</PARAMETER>
<PARAMETER name=\"Sz_total\">1</PARAMETER>
<PARAMETER name=\"MEASURE_AVERAGE[Magnetization]\">Sz</PARAMETER>
<PARAMETER name=\"MEASURE_LOCAL[Local magnetization]\">Sz</PARAMETER>
<PARAMETER name=\"MEASURE_CORRELATIONS[Spin correlations]\">Sz:Sz</PARAMETER>
<PARAMETER name=\"MEASURE_CORRELATIONS[Transverse correlations]\">Splus:Sminus</PARAMETER>
<PARAMETER name=\"MEASURE_STRUCTURE_FACTOR[Structure factor]\">Sz</PARAMETER>
</PARAMETERS></SIMULATION>
")
file(WRITE "${WORK_DIR}/job.in.xml" "<JOB><OUTPUT file=\"job.out.xml\"/>"
  "<TASK status=\"new\"><INPUT file=\"task.in.xml\"/><OUTPUT file=\"task.out.xml\"/></TASK></JOB>\n")

execute_process(COMMAND "${APP}" --write-xml job.in.xml WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "fulldiag exit status: ${result}\nstdout:\n${out}\nstderr:\n${err}")
endif()
file(READ "${WORK_DIR}/task.out.xml" xml)

# sets <var> to the text of <text> from <begin> up to the following <end>
function(between var text begin end)
  string(FIND "${text}" "${begin}" start)
  if(start EQUAL -1)
    message(FATAL_ERROR "no ${begin} in\n${text}")
  endif()
  string(LENGTH "${begin}" length)
  math(EXPR start "${start} + ${length}")
  string(SUBSTRING "${text}" ${start} -1 text)
  string(FIND "${text}" "${end}" stop)
  string(SUBSTRING "${text}" 0 ${stop} text)
  set(${var} "${text}" PARENT_SCOPE)
endfunction()

between(ground "${xml}" "<EIGENSTATE number=\"0\">" "</EIGENSTATE>")

# checks that the ground-state value of <name>, at <index> for a vector
# measurement, lies strictly between <low> and <high>
function(expect name index low high)
  set(text "${ground}")
  if(index STREQUAL "")
    set(key "<SCALAR_AVERAGE name=\"${name}\"><MEAN>")
  else()
    between(text "${text}" "<VECTOR_AVERAGE name=\"${name}\">" "</VECTOR_AVERAGE>")
    set(key "<SCALAR_AVERAGE indexvalue=\"${index}\"><MEAN>")
  endif()
  between(value "${text}" "${key}" "</MEAN>")
  if(NOT (value GREATER ${low} AND value LESS ${high}))
    message(FATAL_ERROR "${name} ${index} = ${value}, expected between ${low} and ${high}\n${ground}")
  endif()
endfunction()

# E = -1/4 - 1/sqrt(2); <Sz_0> = (2 + sqrt(2))/8, <Sz_1> = (2 - sqrt(2))/8;
# <Sz_0 Sz_3> = sqrt(2)/8; <S+_0 S-_3> = -(2 - sqrt(2))/8; S(q = 0) = Sz_total^2/L
expect("Energy" "" -0.957106782 -0.957106780)
expect("Magnetization" "" 0.249999999 0.250000001)
expect("Local magnetization" "( 0 )" 0.426776694 0.426776696)
expect("Local magnetization" "( 1 )" 0.073223304 0.073223306)
expect("Local magnetization" "( 3 )" 0.426776694 0.426776696)
expect("Spin correlations" "( 0 ) -- ( 0 )" 0.249999999 0.250000001)
expect("Spin correlations" "( 0 ) -- ( 3 )" 0.176776694 0.176776696)
expect("Transverse correlations" "( 0 ) -- ( 0 )" 0.926776694 0.926776696)
expect("Transverse correlations" "( 0 ) -- ( 3 )" -0.073223306 -0.073223304)
expect("Structure factor" "( 0 )" 0.249999999 0.250000001)
