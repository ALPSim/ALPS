# n-site <NSITETERM>s in dmrg: Hubbard ladders with a 4-site pair-hopping
# term and a constant. The reference energies are exact diagonalization
# (fulldiag) of the same input, see
# src/alps/model/docs/n_body_terms/benchmarks/alps_checks/t7_dmrg.py.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

file(WRITE "${test_dir}/lattices.xml" [=[<LATTICES>
  <LATTICE name="chain lattice" dimension="1">
    <PARAMETER name="a" default="1"/>
    <BASIS><VECTOR>a</VECTOR></BASIS>
  </LATTICE>
  <UNITCELL name="ladder cell" dimension="1" vertices="2">
    <VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX>
    <VERTEX id="2"><COORDINATE>0.5</COORDINATE></VERTEX>
    <EDGE type="1"><SOURCE vertex="1" offset="0"/><TARGET vertex="2" offset="0"/></EDGE>
    <EDGE type="0"><SOURCE vertex="1" offset="0"/><TARGET vertex="1" offset="1"/></EDGE>
    <EDGE type="0"><SOURCE vertex="2" offset="0"/><TARGET vertex="2" offset="1"/></EDGE>
    <NSITE type="0">
      <SITE vertex="1" offset="0"/>
      <SITE vertex="2" offset="0"/>
      <SITE vertex="2" offset="1"/>
      <SITE vertex="1" offset="1"/>
    </NSITE>
  </UNITCELL>
  <LATTICEGRAPH name="plaquette ladder">
    <FINITELATTICE>
      <LATTICE ref="chain lattice"/>
      <EXTENT dimension="1" size="L"/>
      <BOUNDARY type="open"/>
    </FINITELATTICE>
    <UNITCELL ref="ladder cell"/>
  </LATTICEGRAPH>
</LATTICES>
]=])

file(READ "${xml_dir}/models.xml" models)
string(REPLACE "</MODELS>" [=[
<HAMILTONIAN name="pair hopping ladder">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/><PARAMETER name="K" default="0"/>
  <PARAMETER name="C" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>
  <BONDTERM source="i" target="j">
    -t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))
  </BONDTERM>
  <NSITETERM sites="i j k l">
    -K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))
       + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )
    + C
  </NSITETERM>
</HAMILTONIAN>
<HAMILTONIAN name="spin flip ladder">
  <PARAMETER name="K" default="0"/>
  <BASIS ref="spin"/>
  <BONDTERM source="i" target="j">Sz(i)*Sz(j) + 0.5*(Splus(i)*Sminus(j) + Sminus(i)*Splus(j))</BONDTERM>
  <NSITETERM sites="i j k l">
    K*(Splus(i)*Splus(j)*Splus(k)*Splus(l) + Sminus(i)*Sminus(j)*Sminus(k)*Sminus(l))
  </NSITETERM>
</HAMILTONIAN>
</MODELS>]=] models "${models}")
file(WRITE "${test_dir}/models.xml" "${models}")

file(WRITE "${test_dir}/parameters" [=[
LATTICE_LIBRARY="lattices.xml"
MODEL_LIBRARY="models.xml"
LATTICE="plaquette ladder"
MODEL="pair hopping ladder"
CONSERVED_QUANTUMNUMBERS="Nup,Ndown"
U=8
SWEEPS=4
MAXSTATES=200
NUMBER_EIGENVALUES=1
{ L=3; K=0.7; C=0.3; Nup_total=2; Ndown_total=2; }
{ L=4; K=1.1; C=0.3; Nup_total=3; Ndown_total=2; }
]=])

execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
  WORKING_DIRECTORY "${test_dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)

# fulldiag: -4.766960604357, -6.958608645633
set(expected_1 "-4.76696")
set(expected_2 "-6.95860")
foreach(task 1 2)
  file(READ "${test_dir}/parameters.task${task}.out.xml" result)
  string(REGEX MATCH "<SCALAR_AVERAGE name=\"Energy\"><MEAN>([^<]+)</MEAN>" found "${result}")
  if(NOT found)
    message(FATAL_ERROR "Task ${task}: no energy in the result")
  endif()
  string(FIND "${CMAKE_MATCH_1}" "${expected_${task}}" at)
  if(NOT at EQUAL 0)
    message(FATAL_ERROR "Task ${task}: E0 = ${CMAKE_MATCH_1}, expected ${expected_${task}}...")
  endif()
endforeach()

# A term that changes a quantum number is refused, also when that quantum
# number is not conserved
file(WRITE "${test_dir}/refused" [=[
LATTICE_LIBRARY="lattices.xml"
MODEL_LIBRARY="models.xml"
LATTICE="plaquette ladder"
MODEL="spin flip ladder"
L=3
K=0.5
SWEEPS=2
MAXSTATES=50
{}
]=])
execute_process(COMMAND "${parameter2xml}" refused
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${dmrg}" --write-xml refused.in.xml
  WORKING_DIRECTORY "${test_dir}" RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(status EQUAL 0 OR NOT "${out}${err}" MATCHES "changes the total Sz")
  message(FATAL_ERROR "dmrg did not refuse an n-site term that changes Sz:\n${out}${err}")
endif()
