#  Copyright Lukas Gamper and Synge Todo 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ${cmd} ${binarydir} ${dllexedir})

find_file(input_path ${input}.input ${binarydir} ${sourcedir})
if(NOT input_path)
  find_file(input_path ${input}.ip ${binarydir} ${sourcedir})
endif(NOT input_path)

find_file(output_path ${output}.output ${binarydir} ${sourcedir})
if(NOT output_path)
  find_file(output_path ${output}.op ${binarydir} ${sourcedir})
endif(NOT output_path)

set(ENV{OMP_NUM_THREADS} 1)

if(input_path)
  execute_process(
    COMMAND ${cmd_path}
    RESULT_VARIABLE not_successful
    INPUT_FILE ${input_path}
    OUTPUT_FILE ${cmd}_output
    ERROR_VARIABLE err
    TIMEOUT 600
  )
else(input_path)
  execute_process(
    COMMAND ${cmd_path}
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ${cmd}_output
    ERROR_VARIABLE err
    TIMEOUT 600
  )
endif(input_path)

if(not_successful)
    message(SEND_ERROR "error runing test '${cmd}': ${err};shell output: ${not_successful}!")
endif(not_successful)

if(output_path)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files ${output_path} ${cmd}_output
    RESULT_VARIABLE not_successful
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    TIMEOUT 600
  )
  if(not_successful)
    message(SEND_ERROR "output does not match for '${cmd}': ${err}; ${out}; shell output: ${not_successful}!")
  endif(not_successful)
endif(output_path)

file(REMOVE ${cmd}_output)
