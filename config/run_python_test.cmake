#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

file(WRITE tmp_${cmd}.sh "PYTHONPATH=\$PYTHONPATH:${binarydir}/lib python ${sourcedir}/${cmd}")

execute_process(
    COMMAND sh tmp_${cmd}.sh
    RESULT_VARIABLE not_successful
    ERROR_VARIABLE err
    TIMEOUT 600
)

file(REMOVE tmp_${cmd}.sh)

if(not_successful)
    message(SEND_ERROR "error runing test '${cmd}': ${err}; shell output: ${not_successful}!")
endif(not_successful)
