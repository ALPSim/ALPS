# MEASURE_LOCAL, MEASURE_CORRELATIONS and MEASURE_AVERAGE of operators with
# several terms or a coefficient (issue #304). Sx(i)*Sx(i) and Sx:Sx on the
# spin-1 chain have terms that change Sz, such as S+S+.
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
{
  MODEL="spin"
  local_S=1
  L=5
  Jz=1
  Jxy=0.6
  D=0.4
  h=0.2
  CONSERVED_QUANTUMNUMBERS="Sz"
  Sz_total=1
  MEASURE_LOCAL[sx2]="Sx(i)*Sx(i)"
  MEASURE_AVERAGE[sx2sum]="Sx(i)*Sx(i)"
  MEASURE_CORRELATIONS[sxsx]="Sx:Sx"
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
#             Sx:Sx, rows i = 0..4 over j = 0..4, and Sx(i)*Sx(i) its diagonal,
#             which adds up to 3.144179:
#               0.627272 -0.473753  0.181733 -0.131753  0.094694
#              -0.473753  0.640855 -0.353629  0.155503 -0.131753
#               0.181733 -0.353629  0.607925 -0.353629  0.181733
#              -0.131753  0.155503 -0.353629  0.640855 -0.473753
#               0.094694 -0.131753  0.181733 -0.473753  0.627272
set(task_ex 1)
set(bounds_ex -0.3875:-0.3855 -0.5218:-0.5198 -0.1826:-0.1806 -0.5218:-0.5198 -0.3875:-0.3855)
set(task_exy 1)
set(bounds_exy -0.2911:-0.2891 -0.3704:-0.3684 -0.1692:-0.1672 -0.3704:-0.3684 -0.2911:-0.2891)
set(task_twoSz 1)
set(bounds_twoSz 0.6754:0.6774 -0.0654:-0.0634 0.3870:0.3890 0.3870:0.3890 -0.0654:-0.0634 0.6754:0.6774)
set(task_docc 2)
set(bounds_docc 0.1398:0.1418 0.2295:0.2315 0.2295:0.2315 0.1398:0.1418)
set(task_sx2 3)
set(bounds_sx2 0.6263:0.6283 0.6399:0.6419 0.6069:0.6089 0.6399:0.6419 0.6263:0.6283)
set(task_sxsx 3)
set(bounds_sxsx
  0.6263:0.6283 -0.4748:-0.4728 0.1807:0.1827 -0.1328:-0.1308 0.0937:0.0957
  -0.4748:-0.4728 0.6399:0.6419 -0.3546:-0.3526 0.1545:0.1565 -0.1328:-0.1308
  0.1807:0.1827 -0.3546:-0.3526 0.6069:0.6089 -0.3546:-0.3526 0.1807:0.1827
  -0.1328:-0.1308 0.1545:0.1565 -0.3546:-0.3526 0.6399:0.6419 -0.4748:-0.4728
  0.0937:0.0957 -0.1328:-0.1308 0.1807:0.1827 -0.4748:-0.4728 0.6263:0.6283)

foreach(name ex exy twoSz docc sx2 sxsx)
  file(READ "${test_dir}/parameters.task${task_${name}}.out.xml" result)
  string(REGEX MATCH "<VECTOR_AVERAGE name=\"${name}\"[^>]*>(.*)" found "${result}")
  if(NOT found)
    message(FATAL_ERROR "${name}: no measurement in the result")
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

file(READ "${test_dir}/parameters.task3.out.xml" result)
string(REGEX MATCH "<SCALAR_AVERAGE name=\"sx2sum\"><MEAN>([^<]+)</MEAN>" found "${result}")
if(NOT found)
  message(FATAL_ERROR "sx2sum: no measurement in the result")
endif()
if(CMAKE_MATCH_1 LESS 3.1432 OR CMAKE_MATCH_1 GREATER 3.1452)
  message(FATAL_ERROR "sx2sum = ${CMAKE_MATCH_1}, expected between 3.1432 and 3.1452")
endif()
