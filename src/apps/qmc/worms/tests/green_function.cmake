# The worm has no estimator for the Green's function. A task that asks for it
# must stop with an error that names the parameter, not report zeros.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
file(WRITE "${test_dir}/chain" [=[
LATTICE="chain lattice"
L=4
MODEL="hardcore boson"
t=1
mu=0.2
T=1
THERMALIZATION=100
SWEEPS=500
MEASURE[Green Function]=true
]=])

execute_process(COMMAND "${parameter2xml}" chain
  WORKING_DIRECTORY "${test_dir}" OUTPUT_QUIET RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "parameter2xml failed (${status})")
endif()
execute_process(COMMAND "${worm}" --Tmin 1 chain.in.xml
  WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
  RESULT_VARIABLE status)
if(status EQUAL 0)
  message(FATAL_ERROR "worm accepted MEASURE[Green Function]:\n${out}${err}")
endif()
if(NOT "${out}${err}" MATCHES "MEASURE\\[Green Function\\] is not available in the worm application")
  message(FATAL_ERROR "worm did not say why it stopped:\n${out}${err}")
endif()
