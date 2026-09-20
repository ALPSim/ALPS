#  Copyright Synge Todo 2010-2015.
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

include_guard(GLOBAL)

# Preserve the historical positional API while passing the actual executable
# path to CTest. TARGET_FILE handles suffixes and multi-configuration layouts.
function(add_alps_test name)
  set(program ${name})
  set(input ${name})
  set(output ${name})
  if(ARGC EQUAL 4)
    set(program ${ARGV1})
    set(input ${ARGV2})
    set(output ${ARGV3})
  elseif(ARGC GREATER 1)
    set(input ${ARGV1})
    set(output ${ARGV1})
    if(ARGC EQUAL 3)
      set(output ${ARGV2})
    endif()
  endif()
  add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
    "-Dname=${name}" "-Dcmd_path=$<TARGET_FILE:${program}>"
    "-Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}"
    "-Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}"
    "-Dinput=${input}" "-Doutput=${output}"
    -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/run_test.cmake")
endfunction()

function(add_alps_test_mpi name)
  set(program ${name})
  set(procs 1)
  set(opt "")
  set(input ${name})
  set(output ${name})
  set(args "${ARGN}")
  if(ARGC EQUAL 6)
    list(POP_FRONT args program)
  endif()
  list(LENGTH args count)
  if(count GREATER 0)
    list(POP_FRONT args procs)
  endif()
  if(count GREATER 1)
    list(POP_FRONT args opt)
  endif()
  if(count GREATER 2)
    list(POP_FRONT args input)
    set(output ${input})
  endif()
  if(count GREATER 3)
    list(POP_FRONT args output)
  endif()
  add_test(NAME ${name}-np${procs} COMMAND ${CMAKE_COMMAND}
    "-Dname=${name}-np${procs}" "-Dcmd_path=$<TARGET_FILE:${program}>"
    "-Dopt=${opt}" "-Dmpiexec=${MPIEXEC_EXECUTABLE}"
    "-Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}" "-Dprocs=${procs}"
    "-Dmpiexec_preflags=${MPIEXEC_PREFLAGS}" "-Dmpiexec_postflags=${MPIEXEC_POSTFLAGS}"
    "-Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}"
    "-Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}"
    "-Dinput=${input}" "-Doutput=${output}"
    -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/run_test_mpi.cmake")
endfunction()
