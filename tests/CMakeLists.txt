#  Copyright Matthias Troyer and Synge Todo 2009 - 2012.
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
# Python is optional for native builds, including builds of C++ tests. These
# catalog/CLI tests run when an interpreter and the maintainer dependencies are available.
if(NOT ALPS_CITATION_PYTHON)
  find_package(Python3 3.9 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    set(ALPS_CITATION_PYTHON "${Python3_EXECUTABLE}")
  endif()
endif()
if(ALPS_CITATION_PYTHON)
  execute_process(COMMAND "${ALPS_CITATION_PYTHON}" -c "import yaml, jsonschema"
    RESULT_VARIABLE _alps_citation_dependencies OUTPUT_QUIET ERROR_QUIET)
  if(NOT _alps_citation_dependencies EQUAL 0)
    message(STATUS "Skipping citation wrapper tests: install .github/scripts/citations/requirements.txt in ${ALPS_CITATION_PYTHON}")
    set(ALPS_CITATION_PYTHON "")
  endif()
endif()
if(ALPS_CITATION_PYTHON)
add_test(NAME citation_catalog
  COMMAND "${ALPS_CITATION_PYTHON}" "${PROJECT_SOURCE_DIR}/.github/scripts/generate_citations.py" --check)
set_property(TEST citation_catalog PROPERTY LABELS "utility;citations")
add_test(NAME citation_python
  COMMAND "${ALPS_CITATION_PYTHON}" -m unittest discover
          -s "${PROJECT_SOURCE_DIR}/tests/packaging" -p "*citations.py" -v)
set_property(TEST citation_python PROPERTY LABELS "utility;citations")

foreach(_app dmrg qwl dirloop_sse worm spinmc fulldiag sparsediag hirschfye loop
             dmft hybridization interaction fulldiag_evaluate spinmc_evaluate checksign
             qwl_evaluate worm_evaluate simplemc pevaluate maxent)
  if(TARGET ${_app})
    string(REPLACE "_evaluate" "" _component "${_app}")
    if(_app STREQUAL "loop")
      set(_component "looper")
    elseif(_app MATCHES "^(checksign|simplemc|pevaluate|maxent)$")
      set(_component "framework")
    endif()
    add_test(NAME citation_cli_${_app}
      COMMAND "${ALPS_CITATION_PYTHON}" "${PROJECT_SOURCE_DIR}/tests/utility/check_citation_notice.py"
              "$<TARGET_FILE:${_app}>" "${_component}")
    set_property(TEST citation_cli_${_app} PROPERTY LABELS citations)
    if(ALPS_HAVE_MPI AND MPIEXEC_EXECUTABLE)
      add_test(NAME citation_cli_mpi_${_app}
        COMMAND "${ALPS_CITATION_PYTHON}" "${PROJECT_SOURCE_DIR}/tests/utility/check_citation_notice.py"
                "$<TARGET_FILE:${_app}>" "${_component}"
                "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2" ${MPIEXEC_PREFLAGS})
      set_tests_properties(citation_cli_mpi_${_app} PROPERTIES LABELS citations PROCESSORS 2)
    endif()
  endif()
endforeach()

if(TARGET loop)
  add_test(NAME citation_calculation_loop
    COMMAND "${ALPS_CITATION_PYTHON}" "${PROJECT_SOURCE_DIR}/tests/utility/check_citation_startup.py"
            "$<TARGET_FILE:loop>" "${PROJECT_BINARY_DIR}/src/alps/resources")
  set_property(TEST citation_calculation_loop PROPERTY LABELS citations)
endif()

else()
  message(STATUS "Python citation dependencies unavailable: skipping catalog and CLI wrapper tests")
endif()

if(NOT ALPS_BUILD_LIBS_ONLY)
add_subdirectory(accumulator)
add_subdirectory(alea)
add_subdirectory(fixed_capacity)
add_subdirectory(graph)
add_subdirectory(hdf5)
add_subdirectory(lattice)
add_subdirectory(model)
add_subdirectory(ngs)
add_subdirectory(numeric)
add_subdirectory(osiris)
add_subdirectory(parameter)
add_subdirectory(parapack)
add_subdirectory(parser)
add_subdirectory(random)
add_subdirectory(utility)
endif()
