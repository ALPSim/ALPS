# Two calculations in one process must use their own directories and release
# scratch storage without removing unrelated files or saved solver results.
file(REMOVE_RECURSE "${test_dir}")
if(long_relative_paths)
  # Reproduce a long working directory with short, relative TEMP_DIRECTORY
  # parameters. Only the total path exceeds 255 bytes, not a single component.
  string(REPEAT "nested-" 10 component)
  string(LENGTH "${test_dir}" path_length)
  while(path_length LESS 300)
    string(APPEND test_dir "/${component}")
    string(LENGTH "${test_dir}" path_length)
  endwhile()
endif()
set(first_scratch "${test_dir}/first")
set(second_scratch "${test_dir}/second")
if(long_relative_paths)
  set(first_scratch "first")
  set(second_scratch "second")
endif()
file(MAKE_DIRECTORY "${test_dir}/first" "${test_dir}/second" "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml" "${source_dir}/lib/xml/ALPS.xsl"
  DESTINATION "${test_dir}/xml")
foreach(directory first second)
  file(WRITE "${test_dir}/${directory}/block_ALPS_unrelated" "keep me")
endforeach()
file(WRITE "${test_dir}/parameters" "
LATTICE=\"open chain lattice\"
MODEL=\"spin\"
CONSERVED_QUANTUMNUMBERS=\"N,Sz\"
Sz_total=0
J=1
SWEEPS=2
NUMBER_EIGENVALUES=1
MAXSTATES=20
{ L=8; TEMP_DIRECTORY=\"${first_scratch}\"; }
{ L=8; TEMP_DIRECTORY=\"${second_scratch}\"; }
")
execute_process(COMMAND "${CMAKE_COMMAND}" -E env "ALPS_XML_PATH=${test_dir}/xml"
  "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}" RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "parameter2xml failed: ${output}\n${error}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E env "ALPS_XML_PATH=${test_dir}/xml"
  "${dmrg}" --write-xml parameters.in.xml
  WORKING_DIRECTORY "${test_dir}" RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
file(WRITE "${test_dir}/dmrg.log" "${output}\n${error}")
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dmrg failed: ${output}\n${error}")
endif()
foreach(directory first second)
  file(GLOB remaining RELATIVE "${test_dir}/${directory}" "${test_dir}/${directory}/*")
  if(NOT remaining STREQUAL "block_ALPS_unrelated")
    message(FATAL_ERROR "Unexpected scratch contents in ${directory}: ${remaining}")
  endif()
  file(READ "${test_dir}/${directory}/block_ALPS_unrelated" sentinel)
  if(NOT sentinel STREQUAL "keep me")
    message(FATAL_ERROR "Unrelated file was modified")
  endif()
  string(FIND "${output}" "Creating temp file ${test_dir}/${directory}/" created)
  if(created EQUAL -1)
    message(FATAL_ERROR "No scratch files were created in ${directory}")
  endif()
endforeach()
foreach(task 1 2)
  file(READ "${test_dir}/parameters.task${task}.out.xml" result)
  # Exact eight-site open spin-1/2 Heisenberg ground-state energy is -3.3749325987.
  if(NOT result MATCHES "-3\\.37493")
    message(FATAL_ERROR "Missing expected ground-state energy in task ${task}")
  endif()
endforeach()
