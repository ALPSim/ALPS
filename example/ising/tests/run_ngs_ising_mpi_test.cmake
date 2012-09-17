#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ${cmd} ${binarydir} ${dllexedir})
find_program(mpi_cmd_path ${mpirun} ${binarydir} ${dllexedir})

set(ENV{OMP_NUM_THREADS} 1)

execute_process(
    COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} -T 15 ${sourcedir}/params/param.h5 sim.mpi
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_ising_mpi_output
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_ising_mpi': ${err}; shell output: ${not_successful}!")
endif(not_successful)

execute_process(
    COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} -T 15 --continue ${sourcedir}/params/param.h5 sim.mpi
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ngs_ising_mpi_output_continue
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_ising_mpi': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE dump)
file(REMOVE ngs_ising_mpi_output)
file(REMOVE ngs_ising_mpi_output_continue)

math(EXPR max_rank ${procs}-1)
foreach(rank RANGE ${max_rank})
    file(REMOVE dump.${rank})
endforeach(rank RANGE ${max_rank})
