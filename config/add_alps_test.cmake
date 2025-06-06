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

macro(add_alps_test)
  if(${ARGC} EQUAL 1)
    set(name ${ARGV0})
    set(cmd ${ARGV0})
    set(input ${ARGV0})
    set(output ${ARGV0})
  else(${ARGC} EQUAL 1)
    if(${ARGC} EQUAL 2)
      set(name ${ARGV0})
      set(cmd ${ARGV0})
      set(input ${ARGV1})
      set(output ${ARGV1})
    else(${ARGC} EQUAL 2)
      if(${ARGC} EQUAL 3)
        set(name ${ARGV0})
        set(cmd ${ARGV0})
        set(input ${ARGV1})
        set(output ${ARGV2})
      else(${ARGC} EQUAL 3)
        set(name ${ARGV0})
        set(cmd ${ARGV1})
        set(input ${ARGV2})
        set(output ${ARGV3})
      endif(${ARGC} EQUAL 3)
    endif(${ARGC} EQUAL 2)
  endif(${ARGC} EQUAL 1)
  enable_testing()
  if(MSVC)
    get_target_property(EXE_NAME ${cmd} LOCATION)
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)

  if(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test.cmake)
    set(RUN_TEST ${RUN_TEST_DIR}/run_test.cmake)
  else(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test.cmake)
    if(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test.cmake)
      set(RUN_TEST ${PROJECT_SOURCE_DIR}/config/run_test.cmake)
    else(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test.cmake)
      if(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test.cmake)
        set(RUN_TEST ${ALPS_ROOT_DIR}/share/alps/run_test.cmake)
      else(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test.cmake)
        set(RUN_TEST ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake)
      endif(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test.cmake)
    endif(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test.cmake)
  endif(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test.cmake)
    
  add_test(${name}
    ${CMAKE_COMMAND}
      -Dcmd=${cmd}
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=${input}
      -Doutput=${output}
      -P ${RUN_TEST}
    )
endmacro(add_alps_test)

macro(add_alps_test_mpi)
  if(${ARGC} EQUAL 1)
    set(name ${ARGV0})
    set(cmd ${ARGV0})
    set(procs 1)
    set(opt "")
    set(input ${ARGV0})
    set(output ${ARGV0})
  else(${ARGC} EQUAL 1)
    if(${ARGC} EQUAL 2)
      set(name ${ARGV0})
      set(cmd ${ARGV0})
      set(procs ${ARGV1})
      set(opt "")
      set(input ${ARGV0})
      set(output ${ARGV0})
    else(${ARGC} EQUAL 2)
      if(${ARGC} EQUAL 3)
        set(name ${ARGV0})
        set(mcd ${ARGV0})
        set(procs ${ARGV1})
        set(opt ${ARGV2})
        set(input ${ARGV0})
        set(output ${ARGV0})
      else(${ARGC} EQUAL 3)
        if(${ARGC} EQUAL 4)
          set(name ${ARGV0})
          set(cmd ${ARGV0})
          set(procs ${ARGV1})
          set(opt ${ARGV2})
          set(input ${ARGV3})
          set(output ${ARGV3})
        else(${ARGC} EQUAL 4)
          if(${ARGC} EQUAL 5)
            set(name ${ARGV0})
            set(cmd ${ARGV0})
            set(procs ${ARGV1})
            set(opt ${ARGV2})
            set(input ${ARGV3})
            set(output ${ARGV4})
          else(${ARGC} EQUAL 5)
            set(name ${ARGV0})
            set(cmd ${ARGV1})
            set(procs ${ARGV2})
            set(opt ${ARGV3})
            set(input ${ARGV4})
            set(output ${ARGV5})
          endif(${ARGC} EQUAL 5)
        endif(${ARGC} EQUAL 4)
      endif(${ARGC} EQUAL 3)
    endif(${ARGC} EQUAL 2)
  endif(${ARGC} EQUAL 1)
  enable_testing()
  if(MSVC)
    get_target_property(EXE_NAME ${cmd} LOCATION)
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
  
  if(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test_mpi.cmake)
    set(RUN_TEST ${RUN_TEST_DIR}/run_test_mpi.cmake)
  else(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test_mpi.cmake)
    if(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake)
      set(RUN_TEST ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake)
    else(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake)
      if(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test_mpi.cmake)
        set(RUN_TEST ${ALPS_ROOT_DIR}/share/alps/run_test_mpi.cmake)
      else(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test_mpi.cmake)
        set(RUN_TEST ${CMAKE_INSTALL_PREFIX}/share/alps/run_test_mpi.cmake)
      endif(EXISTS ${ALPS_ROOT_DIR}/share/alps/run_test_mpi.cmake)
    endif(EXISTS ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake)
  endif(RUN_TEST_DIR AND EXISTS ${RUN_TEST_DIR}/run_test_mpi.cmake)
    
  add_test(${name}-np${procs}
    ${CMAKE_COMMAND}
      -Dcmd=${cmd}
      -Dopt=${opt}
      -Dmpiexec=${MPIEXEC}
      -Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}
      -Dprocs=${procs}
      -Dmpiexec_preflags=${MPIEXEC_PREFLAGS}
      -Dmpiexec_postflags=${MPIEXEC_POSTFLAGS}
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=${input}
      -Doutput=${output}
      -P ${RUN_TEST}
    )
endmacro(add_alps_test_mpi)
