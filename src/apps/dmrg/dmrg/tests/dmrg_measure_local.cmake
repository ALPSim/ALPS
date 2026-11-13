# MEASURE_LOCAL of operators with several terms or a coefficient (issue #304).
# The reference values are sparsediag of the same input; each value has to
# lie within 0.001 of it, given below as lower:upper bounds.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

file(WRITE "${test_dir}/parameters" [=[
LATTICE="open chain lattice"
NUMBER_EIGENVALUES=1
SWEEPS=6
MAXSTATES=100
{
  MODEL="spin"
  local_S=1/2
  L=6
  Jz=0.7
  Jxy=1
  h=0.3
  CONSERVED_QUANTUMNUMBERS="Sz"
  Sz_total=1
  MEASURE_LOCAL[ex]="exchange"
  MEASURE_LOCAL[exy]="exchange_xy"
  MEASURE_LOCAL[twoSz]="2*Sz"
}
{
  MODEL="boson Hubbard"
  Nmax=2
  L=4
  t=1
  U=1.5
  CONSERVED_QUANTUMNUMBERS="N"
  N_total=4
  MEASURE_LOCAL[docc]="double_occupancy"
}
]=])

execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
  WORKING_DIRECTORY "${test_dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)

# sparsediag: exchange -0.386488 -0.520837 -0.181568 -0.520837 -0.386488
#             exchange_xy -0.290123 -0.369372 -0.168166 -0.369372 -0.290123
#             2*Sz 0.676394 -0.064401 0.388007 0.388007 -0.064401 0.676394
#             double_occupancy 0.140755 0.230531 0.230531 0.140755
set(task_ex 1)
set(bounds_ex -0.3875:-0.3855 -0.5218:-0.5198 -0.1826:-0.1806 -0.5218:-0.5198 -0.3875:-0.3855)
set(task_exy 1)
set(bounds_exy -0.2911:-0.2891 -0.3704:-0.3684 -0.1692:-0.1672 -0.3704:-0.3684 -0.2911:-0.2891)
set(task_twoSz 1)
set(bounds_twoSz 0.6754:0.6774 -0.0654:-0.0634 0.3870:0.3890 0.3870:0.3890 -0.0654:-0.0634 0.6754:0.6774)
set(task_docc 2)
set(bounds_docc 0.1398:0.1418 0.2295:0.2315 0.2295:0.2315 0.1398:0.1418)

foreach(name ex exy twoSz docc)
  file(READ "${test_dir}/parameters.task${task_${name}}.out.xml" result)
  string(REGEX MATCH "<VECTOR_AVERAGE name=\"${name}\"[^>]*>(.*)" found "${result}")
  if(NOT found)
    message(FATAL_ERROR "${name}: no local measurement in the result")
  endif()
  string(REGEX REPLACE "</VECTOR_AVERAGE>.*" "" block "${CMAKE_MATCH_1}")
  string(REGEX MATCHALL "<MEAN>[^<]+</MEAN>" means "${block}")
  list(LENGTH means n)
  list(LENGTH bounds_${name} n_expected)
  if(NOT n EQUAL n_expected)
    message(FATAL_ERROR "${name}: ${n} values, expected ${n_expected}")
  endif()
  set(i 0)
  foreach(bound IN LISTS bounds_${name})
    list(GET means ${i} mean)
    string(REGEX REPLACE "</?MEAN>" "" value "${mean}")
    string(REPLACE ":" ";" range "${bound}")
    list(GET range 0 lower)
    list(GET range 1 upper)
    if(value LESS lower OR value GREATER upper)
      message(FATAL_ERROR "${name}[${i}] = ${value}, expected between ${lower} and ${upper}")
    endif()
    math(EXPR i "${i} + 1")
  endforeach()
endforeach()
