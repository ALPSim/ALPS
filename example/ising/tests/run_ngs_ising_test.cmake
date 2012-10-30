#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_program(cmd_path ${cmd} ${binarydir} ${bindir})

set(ENV{OMP_NUM_THREADS} 1)

execute_process(
    COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} -T 15 ${sourcedir}/params/param.h5 ${cmd}.result
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ${cmd}_output
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test '${cmd}': ${err}; shell output: ${not_successful}!")
endif(not_successful)

execute_process(
    COMMAND ${mpiexec} ${mpiexec_numproc_flag} ${procs} ${cmd_path} -T 15 --continue ${sourcedir}/params/param.h5 ${cmd}.result
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ${cmd}_output_continue
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test '${cmd}': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE ${cmd}.result)
file(REMOVE ${cmd}_output)
file(REMOVE ${cmd}_output_continue)

file(GLOB checkpoint_files checkpoint*)
foreach (name ${checkpoint_files})
    file(REMOVE ${name})
endforeach(name)
