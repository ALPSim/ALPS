# Two small tasks exercise the solver with separate scratch directories.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/scratch1" "${test_dir}/scratch2" "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml"
  "${source_dir}/src/alps/resources/ALPS.xsl" DESTINATION "${test_dir}/xml")
set(ENV{ALPS_XML_PATH} "${test_dir}/xml")

file(WRITE "${test_dir}/parameters" [=[
LATTICE="open chain lattice"
MODEL="spin"
CONSERVED_QUANTUMNUMBERS="N,Sz"
Sz_total=0
L=8
J=1
SWEEPS=2
NUMBER_EIGENVALUES=1
MAXSTATES=20
]=])
foreach(task 1 2)
  file(WRITE "${test_dir}/scratch${task}/block_ALPS_unrelated" "keep me")
  file(APPEND "${test_dir}/parameters"
    "{ TEMP_DIRECTORY=\"${test_dir}/scratch${task}\"; }\n")
endforeach()

execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
  WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)

foreach(task 1 2)
  if(NOT EXISTS "${test_dir}/parameters.task${task}.out.xml")
    message(FATAL_ERROR "Task ${task}: result file missing")
  endif()
  file(GLOB remaining RELATIVE "${test_dir}/scratch${task}"
    "${test_dir}/scratch${task}/*")
  if(NOT remaining STREQUAL "block_ALPS_unrelated")
    message(FATAL_ERROR "Task ${task}: unexpected scratch contents: ${remaining}")
  endif()
  file(READ "${test_dir}/scratch${task}/block_ALPS_unrelated" sentinel)
  if(NOT sentinel STREQUAL "keep me")
    message(FATAL_ERROR "Task ${task}: unrelated file was modified")
  endif()
endforeach()
