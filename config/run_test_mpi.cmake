#  Copyright Lukas Gamper and Synge Todo 2009 - 2011.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

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
