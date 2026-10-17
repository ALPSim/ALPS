# SAVE_MPS=1 writes the final ground state as /mps; dmrg_mps_check reads it
# back. The spin chain closes its final sweep at L/2, the Hubbard chain, with
# two targets and correlations (stored products), at L-2.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml"
  "${xml_dir}/ALPS.xsl" DESTINATION "${test_dir}/xml")
set(ENV{ALPS_XML_PATH} "${test_dir}/xml")

file(WRITE "${test_dir}/parameters" [=[
LATTICE="open chain lattice"
SWEEPS=4
MAXSTATES=64
SAVE_MPS=1
{ MODEL="spin"; CONSERVED_QUANTUMNUMBERS="Sz"; L=8; J=1; Sz_total=1; }
{ MODEL="fermion Hubbard"; CONSERVED_QUANTUMNUMBERS="Nup,Ndown"; L=8; t=1; U=4; Nup_total=3; Ndown_total=3;
  NUMBER_EIGENVALUES=2; MEASURE_CORRELATIONS[nn]="n:n"; }
]=])

execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
  WORKING_DIRECTORY "${test_dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)

foreach(check "1;Sz;1;4" "2;Nup;3;6")
  list(GET check 0 task)
  list(GET check 1 charge)
  list(GET check 2 total)
  list(GET check 3 center)
  execute_process(COMMAND "${checker}" "parameters.task${task}.out.h5" "${charge}" "${total}"
    WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Task ${task}: the saved MPS failed its check")
  endif()
  if(NOT out MATCHES "center ${center}\n")
    message(FATAL_ERROR "Task ${task}: expected the center at site ${center}: ${out}")
  endif()
endforeach()
