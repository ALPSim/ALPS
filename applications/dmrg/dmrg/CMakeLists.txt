#  Copyright Matthias Troyer 2009 - 2010.
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

if(LAPACK_FOUND)
if(ALPS_LLVM_WORKAROUND)
  message(STATUS "dmrg will not be built since llvm workaround is enabled")
else(ALPS_LLVM_WORKAROUND)
  add_definitions(${LAPACK_DEFINITIONS} -DBIND_FORTRAN_LOWERCASE_UNDERSCORE)
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmrg/dmrg)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(dmrg dmrg.C factory.C)
  target_link_libraries(dmrg alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
if(ALPS_PYTHON_WHEEL)
  install(TARGETS dmrg RUNTIME DESTINATION pyalps/bin COMPONENT applications)
set_target_properties(dmrg  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
else()
  install(TARGETS dmrg RUNTIME DESTINATION bin COMPONENT applications)
endif()
endif(ALPS_LLVM_WORKAROUND)
else(LAPACK_FOUND)
  message(STATUS "dmrg will not be built since lapack library is not found")
endif(LAPACK_FOUND)
