# Site and bond terms that change a quantum number of the site basis are
# refused (issue #296): dmtk dropped them without a message.
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

function(expect_refusal name message)
  run_dmrg(${name})
  if(status EQUAL 0 OR NOT output MATCHES "${message}")
    message(FATAL_ERROR "${name}: dmrg did not refuse with \"${message}\":\n${output}")
  endif()
endfunction()

# Transverse field -Gamma*Sx(i), as in the issue
file(WRITE "${test_dir}/field" [=[
LATTICE="open chain lattice"
MODEL="spin"
local_S=1/2
L=8
SWEEPS=2
MAXSTATES=20
{ Jz=1; Jxy=1; Gamma=0.5; }
]=])
expect_refusal(field "A <SITETERM> changes Sz")

# Also with Sz conserved
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
expect_refusal(conserved "A <SITETERM> changes Sz")

# A bond term that changes Sz by 2
file(WRITE "${test_dir}/xising" [=[
MODEL_LIBRARY="models.xml"
LATTICE="open chain lattice"
MODEL="x Ising"
local_S=1/2
L=8
SWEEPS=2
MAXSTATES=20
{}
]=])
expect_refusal(xising "A <BONDTERM> changes Sz")

# A pair field changes Nup and Ndown
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
expect_refusal(pairs "A <SITETERM> changes N")

# A field along z conserves Sz and is kept: sparsediag -3.374932598688
file(WRITE "${test_dir}/zfield" [=[
LATTICE="open chain lattice"
MODEL="spin"
local_S=1/2
L=8
SWEEPS=4
MAXSTATES=60
NUMBER_EIGENVALUES=1
{ Jz=1; Jxy=1; h=0.3; }
]=])
run_dmrg(zfield)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dmrg failed on a field along z:\n${output}")
endif()
file(READ "${test_dir}/zfield.task1.out.xml" result)
string(REGEX MATCH "<SCALAR_AVERAGE name=\"Energy\"><MEAN>([^<]+)</MEAN>" found "${result}")
string(FIND "${CMAKE_MATCH_1}" "-3.37493" at)
if(NOT found OR NOT at EQUAL 0)
  message(FATAL_ERROR "Field along z: E0 = ${CMAKE_MATCH_1}, expected -3.37493...")
endif()
