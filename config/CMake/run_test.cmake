find_program(cmd_path ${cmd} ${binarydir} ${dllexedir})
find_file(input_path ${input} ${binarydir} ${sourcedir})
find_file(output_path ${output} ${binarydir} ${sourcedir})

if(EXISTS ${input_path})
  execute_process(
	COMMAND ${cmd_path}
	RESULT_VARIABLE not_successful
	INPUT_FILE ${input_path}
	OUTPUT_FILE ${cmd}_output
	ERROR_VARIABLE err
	TIMEOUT ${timeout}
  )
else(EXISTS ${input_path})
  execute_process(
	COMMAND ${cmd_path}
	RESULT_VARIABLE not_successful
	OUTPUT_FILE ${cmd}_output
	ERROR_VARIABLE err
	TIMEOUT ${timeout}
  )
endif(EXISTS ${input_path})

if(not_successful)
	message(SEND_ERROR "error runing test '${cmd}': ${err};shell output: ${not_successful}!")
endif(not_successful)

# file(WRITE ${cmd}_output ${out})

if(EXISTS ${output_path})
  execute_process(
	COMMAND ${CMAKE_COMMAND} -E compare_files ${output_path} ${cmd}_output
	RESULT_VARIABLE not_successful
	OUTPUT_VARIABLE out
	ERROR_VARIABLE err
	TIMEOUT ${timeout}
  )
  if(not_successful)
  	message(SEND_ERROR "output does not match for '${cmd}': ${err}; ${out}; shell output: ${not_successful}!")
  endif(not_successful)
endif(EXISTS ${output_path})

#file(REMOVE ${cmd}_output)
