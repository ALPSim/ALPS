#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
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

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

#add_executable(loadobs loadobs.cpp)
#add_dependencies(loadobs alps)
#target_link_libraries(loadobs alps)
#add_alps_test(loadobs)

enable_testing()
if (ALPS_BUILD_PYTHON AND BUILD_SHARED_LIBS)
#    FOREACH (name pyioarchive pyhdf5io numpylarge pyparams hlist_test mcdata pyhdf5 mcanalyze, accumulators)
    FOREACH (name pyioarchive pyhdf5io_test pyparams_test hlist_test mcdata_test pyhdf5_test mcanalyze)
        add_test(python_${name}
            ${CMAKE_COMMAND}
                -Dpython_interpreter=${PYTHON_INTERPRETER}
                -Dcmd=${name}.py
                -Dinput=${name}
                -Doutput=${name}
                -Dpythonpath=${PROJECT_BINARY_DIR}/lib/pyalps:${PROJECT_SOURCE_DIR}/lib
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dcmddir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/run_python_test.cmake
        )
        set_property(TEST python_${name} PROPERTY LABELS pyalps)
    ENDFOREACH(name)
ENDIF(ALPS_BUILD_PYTHON AND BUILD_SHARED_LIBS)
