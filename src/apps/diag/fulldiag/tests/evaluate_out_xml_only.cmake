# fulldiag_evaluate reads a task's .out.xml file (#433). Given the .out.h5
# file written next to it, it stops with a usage error naming .out.xml and
# writes nothing; the .out.xml file of the same run still evaluates. The usage
# text names the .out.xml file.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
file(WRITE "${test_dir}/fd" [=[
LATTICE="chain lattice"
L=4
MODEL="spin"
J=1
CONSERVED_QUANTUMNUMBERS="Sz"
]=])

function(run name)
  execute_process(COMMAND ${ARGN}
    WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
    RESULT_VARIABLE status)
  set(status "${status}" PARENT_SCOPE)
  set(output "${out}${err}" PARENT_SCOPE)
endfunction()

run(fulldiag_evaluate "${fulldiag_evaluate}" --help)
if(NOT status EQUAL 0 OR NOT output MATCHES "\\.out\\.xml")
  message(FATAL_ERROR "fulldiag_evaluate --help does not name the .out.xml file "
                      "(${status}):\n${output}")
endif()

run(parameter2xml "${parameter2xml}" fd)
run(fulldiag "${fulldiag}" fd.in.xml)
if(NOT status EQUAL 0 OR NOT EXISTS "${test_dir}/fd.task1.out.h5")
  message(FATAL_ERROR "fulldiag failed (${status}):\n${output}")
endif()

run(fulldiag_evaluate "${fulldiag_evaluate}" --T_MIN 0.5 --T_MAX 1.0 --DELTA_T 0.5 fd.task1.out.h5)
if(NOT status EQUAL 1 OR NOT output MATCHES "\\.out\\.xml")
  message(FATAL_ERROR "fulldiag_evaluate on the .out.h5 file: expected a usage error "
                      "naming .out.xml, got status ${status}:\n${output}")
endif()
file(GLOB plots "${test_dir}/*.plot.*")
if(plots)
  message(FATAL_ERROR "fulldiag_evaluate on the .out.h5 file wrote ${plots}")
endif()

run(fulldiag_evaluate "${fulldiag_evaluate}" --T_MIN 0.5 --T_MAX 1.0 --DELTA_T 0.5 fd.task1.out.xml)
if(NOT status EQUAL 0 OR NOT EXISTS "${test_dir}/fd.task1.plot.energy.xml")
  message(FATAL_ERROR "fulldiag_evaluate on the .out.xml file failed (${status}):\n${output}")
endif()
