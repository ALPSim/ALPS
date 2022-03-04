#  Copyright Lukas Gamper and Synge Todo 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ${cmd} ${binarydir} ${dllexedir} NO_SYSTEM_ENVIRONMENT_PATH)
find_program(mpi_cmd_path ${mpirun} ${binarydir} ${dllexedir} NO_SYSTEM_ENVIRONMENT_PATH)

find_file(input_path ${input}.input-${procs} ${binarydir} ${sourcedir})
if(NOT input_path)
  find_file(input_path ${input}.ip-${procs} ${binarydir} ${sourcedir})
endif(NOT input_path)
if(NOT input_path)
  find_file(input_path ${input}.input ${binarydir} ${sourcedir})
endif(NOT input_path)
if(NOT input_path)
  find_file(input_path ${input}.ip ${binarydir} ${sourcedir})
endif(NOT input_path)

find_file(output_path ${output}.output-${procs} ${binarydir} ${sourcedir})
if(NOT output_path)
  find_file(output_path ${output}.op-${procs} ${binarydir} ${sourcedir})
endif(NOT output_path)

set(ENV{OMP_NUM_THREADS} 1)

if(input_path)
  execute_process(
	COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} ${opt}
	RESULT_VARIABLE not_successful
	INPUT_FILE ${input_path}
	OUTPUT_FILE ${cmd}_output_${procs}
	ERROR_VARIABLE err
	TIMEOUT 600
  )
else(input_path)
  execute_process(
	COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} ${opt}
	RESULT_VARIABLE not_successful
	OUTPUT_FILE ${cmd}_output_${procs}
	ERROR_VARIABLE err
	TIMEOUT 600
  )
endif(input_path)

if(not_successful)
	message(SEND_ERROR "error runing test '${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd} ${opt}': ${err};shell output: ${not_successful}!")
endif(not_successful)

if(output_path)
  execute_process(
	COMMAND ${CMAKE_COMMAND} -E compare_files ${output_path} ${cmd}_output_${procs}
	RESULT_VARIABLE not_successful
	OUTPUT_VARIABLE out
	ERROR_VARIABLE err
	TIMEOUT 600
  )
  if(not_successful)
  	message(SEND_ERROR "output does not match for '${cmd} ${opt}': ${err}; ${out}; shell output: ${not_successful}!")
  endif(not_successful)
endif(output_path)

file(REMOVE ${cmd}_output_${procs})
