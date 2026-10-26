# MEASURE_ENERGY=false without custom measurements writes no per-sector
# results. The file must still be readable by fulldiag_evaluate and convert2xml.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
file(WRITE "${test_dir}/noenergy" [=[
LATTICE="chain lattice"
L=6
MODEL="spin"
J=1
CONSERVED_QUANTUMNUMBERS="Sz"
T=0.5
h=0.3
MEASURE_ENERGY=false
]=])

function(run name)
  execute_process(COMMAND ${ARGN}
    WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${name} failed (${status}):\n${out}${err}")
  endif()
  set(output "${out}${err}" PARENT_SCOPE)
endfunction()

run(parameter2xml "${parameter2xml}" noenergy)
run(fulldiag "${fulldiag}" noenergy.in.xml)
run(fulldiag_evaluate "${fulldiag_evaluate}" --T_MIN 0.1 --T_MAX 1 --DELTA_T 0.3 noenergy.task1.out.xml)
if(NOT EXISTS "${test_dir}/noenergy.task1.plot.energy.xml")
  message(FATAL_ERROR "fulldiag_evaluate wrote no energy plot:\n${output}")
endif()
run(convert2xml "${convert2xml}" noenergy.task1.out.xml)
file(READ "${test_dir}/noenergy.task1.out.xml" task)
if(NOT task MATCHES "<EIGENVALUES")
  message(FATAL_ERROR "convert2xml wrote no eigenvalues:\n${task}")
endif()
