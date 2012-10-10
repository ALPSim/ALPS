#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ngs_ising_single ${binarydir} ${dllexedir})

execute_process(
    COMMAND ${cmd_path} -T 15 ${currentdir}/params/param.h5 sim.sgl
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_ising_sgl_output
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_ising_single': ${err}; shell output: ${not_successful}!")
endif(not_successful)

execute_process(
    COMMAND ${cmd_path} -T 15 --continue ${currentdir}/params/param.h5 sim.sgl
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_ising_sgl_output_continue
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_ising_single': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE dump)
file(REMOVE sim.sgl)
file(REMOVE ngs_ising_sgl_output)
file(REMOVE ngs_ising_sgl_output_continue)
