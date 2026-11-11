# Site and bond terms that change a quantum number of the site basis (issue
# #296). The reference energies are sparsediag of the same input.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

file(READ "${xml_dir}/models.xml" models)
string(REPLACE "</MODELS>" [=[
<HAMILTONIAN name="x Ising">
  <PARAMETER name="J" default="1"/>
  <BASIS ref="spin"/>
  <BONDTERM source="i" target="j">-J*Sx(i)*Sx(j)</BONDTERM>
</HAMILTONIAN>
<HAMILTONIAN name="pair field">
  <PARAMETER name="t" default="1"/><PARAMETER name="D" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">-D*(cdag_up(i)*cdag_down(i) + c_down(i)*c_up(i))</SITETERM>
  <BONDTERM source="i" target="j">
    -t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))
  </BONDTERM>
</HAMILTONIAN>
</MODELS>]=] models "${models}")
file(WRITE "${test_dir}/models.xml" "${models}")

function(run_dmrg name)
  execute_process(COMMAND "${parameter2xml}" ${name}
    WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
  execute_process(COMMAND "${dmrg}" --write-xml ${name}.in.xml
    WORKING_DIRECTORY "${test_dir}" RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(status "${status}" PARENT_SCOPE)
  set(output "${out}${err}" PARENT_SCOPE)
endfunction()

function(check_energy name task expected)
  file(READ "${test_dir}/${name}.task${task}.out.xml" result)
  string(REGEX MATCH "<SCALAR_AVERAGE name=\"Energy\"><MEAN>([^<]+)</MEAN>" found "${result}")
  string(FIND "${CMAKE_MATCH_1}" "${expected}" at)
  if(NOT found OR NOT at EQUAL 0)
    message(FATAL_ERROR "${name} task ${task}: E0 = ${CMAKE_MATCH_1}, expected ${expected}...")
  endif()
endfunction()

# Transverse field -Gamma*Sx(i): sparsediag -2.459487861948, -3.482240487823
file(WRITE "${test_dir}/field" [=[
LATTICE="open chain lattice"
MODEL="spin"
local_S=1/2
L=8
SWEEPS=4
MAXSTATES=60
NUMBER_EIGENVALUES=1
{ Jz=1; Jxy=0; Gamma=0.5; }
{ Jz=1; Jxy=1; Gamma=0.5; }
]=])
run_dmrg(field)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dmrg failed on a transverse field:\n${output}")
endif()
check_energy(field 1 "-2.45948")
check_energy(field 2 "-3.48224")

# A bond term that changes Sz by 2: seven Ising bonds along x, E0 = -7/4
file(WRITE "${test_dir}/xising" [=[
MODEL_LIBRARY="models.xml"
LATTICE="open chain lattice"
MODEL="x Ising"
local_S=1/2
L=8
SWEEPS=4
MAXSTATES=60
NUMBER_EIGENVALUES=1
{}
]=])
run_dmrg(xising)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dmrg failed on an Sx*Sx bond term:\n${output}")
endif()
check_energy(xising 1 "-1.75")

# The field changes Sz, which CONSERVED_QUANTUMNUMBERS fixes
file(WRITE "${test_dir}/conserved" [=[
LATTICE="open chain lattice"
MODEL="spin"
local_S=1/2
L=8
SWEEPS=2
MAXSTATES=20
CONSERVED_QUANTUMNUMBERS="Sz"
Sz_total=0
{ Jz=1; Jxy=1; Gamma=0.5; }
]=])
run_dmrg(conserved)
if(status EQUAL 0 OR NOT output MATCHES "changes Sz, which CONSERVED_QUANTUMNUMBERS fixes")
  message(FATAL_ERROR "dmrg did not refuse a field that changes the conserved Sz:\n${output}")
endif()

# A pair field changes Nup and Ndown, which are fermionic
file(WRITE "${test_dir}/pairs" [=[
MODEL_LIBRARY="models.xml"
LATTICE="open chain lattice"
MODEL="pair field"
L=4
D=0.5
SWEEPS=2
MAXSTATES=20
{}
]=])
run_dmrg(pairs)
if(status EQUAL 0 OR NOT output MATCHES "changes the fermionic quantum number")
  message(FATAL_ERROR "dmrg did not refuse a term that changes a fermionic quantum number:\n${output}")
endif()
