#  Copyright Synge Todo 2010.
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

add_executable(exchange main.C ising.C loop.C)
target_link_libraries(exchange PRIVATE ALPS::alps)

alps_add_test(exchange)

# Under MPI, check the parallel exchange workers and the parallel Ising worker
# against the exact Ising chain (see check_mpi.py).
if(ALPS_BUILD_TESTING AND ALPS_ENABLE_MPI)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND AND MPIEXEC_EXECUTABLE)
    add_test(NAME exchange_mpi
      COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/check_mpi.py" "${CMAKE_CURRENT_SOURCE_DIR}"
        "$<TARGET_FILE:exchange>" "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG})
    set_tests_properties(exchange_mpi PROPERTIES PROCESSORS 4 TIMEOUT 300 LABELS "mpi"
      ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${ALPS_DATA_DIR}/xml;OMP_NUM_THREADS=set:1")
  endif()
endif()
