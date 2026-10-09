#  Copyright Synge Todo 2012.
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

FOREACH(name vmusage bitops citations)
  add_executable(${name} ${name}.cpp)
  add_dependencies(${name} alps)
  target_link_libraries(${name} alps)
  add_alps_test(${name})
  set_property(TEST ${name} PROPERTY LABELS utility)
ENDFOREACH(name)

if(ALPS_CITATION_PYTHON)
add_executable(citation_runtime citation_runtime.cpp)
target_link_libraries(citation_runtime alps)
add_test(NAME citation_runtime_serial
  COMMAND "${ALPS_CITATION_PYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/check_citation_runtime.py"
          "$<TARGET_FILE:citation_runtime>")
set_property(TEST citation_runtime_serial PROPERTY LABELS citations)
if(ALPS_HAVE_MPI AND MPIEXEC_EXECUTABLE)
  add_test(NAME citation_runtime_mpi
    COMMAND "${ALPS_CITATION_PYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/check_citation_runtime.py"
            "$<TARGET_FILE:citation_runtime>"
            "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2" ${MPIEXEC_PREFLAGS})
  set_tests_properties(citation_runtime_mpi PROPERTIES LABELS citations PROCESSORS 2)
endif()

endif()

# filelist unit test
include(${CMAKE_CURRENT_SOURCE_DIR}/dmrg_filelist_tests.cmake)
