#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(mpi_cmd_path ${mpirun} ${binarydir} ${dllexedir})

file(WRITE run_ngs_ising_python_native_tmp.sh "PYTHONPATH=\$PYTHONPATH:${project_binarydir}/lib/pyalps:${project_binarydir}/src/boost:${project_sourcedir}/lib ${mpi_cmd_path} -n 2 ${python_interpreter} ${sourcedir}/main_native.py 15 f sim.mpi.h5")

set(ENV{OMP_NUM_THREADS} 1)

execute_process(
    COMMAND sh run_ngs_ising_python_native_tmp.sh
    RESULT_VARIABLE not_successful
    OUTPUT_FILE run_ngs_ising_python_native_output
    ERROR_VARIABLE err
    TIMEOUT 600
)

file(REMOVE run_ngs_ising_python_native_tmp.sh)
file(REMOVE sim.mpi.h5)

if(not_successful)
    message(SEND_ERROR "error runing test 'ngs_ising_python_native_mpi': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE dump)
file(REMOVE run_ngs_ising_python_native_output)

math(EXPR max_rank ${procs}-1)
foreach(rank RANGE ${max_rank})
    file(REMOVE dump.${rank})
endforeach(rank RANGE ${max_rank})