#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

file(WRITE ${cmd}_tmp.sh "PYTHONPATH=\$PYTHONPATH:${binarydir}/lib/pyalps:${binarydir}/src/boost:${sourcedir}/lib:${bindir}:${currentdir} ${runexec} ${python_interpreter} ${currentdir}/${cmd}.py -T 15 ${cmd}.h5")

execute_process(
    COMMAND sh ${cmd}_tmp.sh
    RESULT_VARIABLE not_successful
    OUTPUT_FILE ${cmd}_output
    ERROR_VARIABLE err
    TIMEOUT 600
)

if(not_successful)
    message(SEND_ERROR "error runing test '${cmd}': ${err}; shell output: ${not_successful}!")
endif(not_successful)

file(REMOVE ${cmd}_tmp.sh)
file(REMOVE ${cmd}.h5)
file(REMOVE ${cmd}_output)
