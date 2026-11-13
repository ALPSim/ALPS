# A <GRAPH> lattice has no periodic lattice directions to count winding
# numbers along. The worm must run anyway, say that
# it skips the winding numbers and the stiffness, and produce the other results.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
file(WRITE "${test_dir}/lattices.xml" [=[
<LATTICES>
<GRAPH name="triangle" vertices="3">
  <EDGE type="0" source="1" target="2"/>
  <EDGE type="0" source="2" target="3"/>
  <EDGE type="0" source="3" target="1"/>
</GRAPH>
</LATTICES>
]=])
file(WRITE "${test_dir}/triangle" [=[
LATTICE_LIBRARY="lattices.xml"
LATTICE="triangle"
MODEL="boson Hubbard"
t=1
U=4
mu=1
Nmax=2
T=0.5
THERMALIZATION=100
SWEEPS=500
SEED=7
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

run(parameter2xml "${parameter2xml}" triangle)
run(worm "${worm}" --Tmin 1 triangle.in.xml)
if(NOT output MATCHES "no periodic lattice directions: winding numbers and the stiffness are not measured")
  message(FATAL_ERROR "worm did not say that it skips the winding numbers:\n${output}")
endif()
run(convert2xml "${convert2xml}" triangle.task1.out.xml)
file(READ "${test_dir}/triangle.task1.out.xml" task)
if(NOT task MATCHES "name=\"Energy\"")
  message(FATAL_ERROR "the result has no energy:\n${task}")
endif()
if(task MATCHES "Winding|Stiffness")
  message(FATAL_ERROR "the result has winding numbers or a stiffness:\n${task}")
endif()
