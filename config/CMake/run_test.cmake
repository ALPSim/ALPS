
find_program(cmd_path ${cmd})

execute_process(
	COMMAND ${cmd_path}
	RESULT_VARIABLE not_successful
	INPUT_FILE ${input}
	OUTPUT_VARIABLE out
	ERROR_VARIABLE err
	TIMEOUT ${timeout}
)

file(WRITE ${cmd}_output ${out})

if(not_successful)
	message(SEND_ERROR "error runing test '${cmd}': ${err};shell output: ${not_successful}!")
endif(not_successful)

execute_process(
	COMMAND ${CMAKE_COMMAND} -E compare_files ${output} ${cmd}_output
	RESULT_VARIABLE not_successful
	OUTPUT_QUIET
	ERROR_VARIABLE err
	TIMEOUT ${timeout}
)

file(REMOVE ${cmd}_output)

if(not_successful)
	message(SEND_ERROR "output does not match for '${cmd}': ${err}; shell output: ${not_successful}!")
endif(not_successful)
