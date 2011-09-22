#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ngs_mpi ${binarydir} ${dllexedir})

execute_process(
    COMMAND mpirun -n 4 ${cmd_path} -T 15 --mpi ${currentdir}/param.h5 sim.pscan
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_mpi_output_1
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_mpi': ${err}; shell output: ${not_successful}!")
endif(not_successful)

execute_process(
    COMMAND mpirun -n 4 ${cmd_path} -T 15 --continue --mpi ${currentdir}/param.h5 sim.pscan
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_mpi_output_2
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_mpi': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE ngs_mpi_output_1)
file(REMOVE ngs_mpi_output_2)
