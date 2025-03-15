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
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_library(fulldiag_impl STATIC factory.C)
  target_link_libraries(fulldiag_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

  add_executable(fulldiag fulldiag.C)
  add_executable(fulldiag_evaluate fulldiag_evaluate.C)
  target_link_libraries(fulldiag fulldiag_impl)
  target_link_libraries(fulldiag_evaluate fulldiag_impl)
  if(ALPS_PYTHON_WHEEL)
  if(APPLE)
    set_target_properties(fulldiag PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
    set_target_properties(fulldiag_evaluate PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
  else(APPLE)
    set_target_properties(fulldiag PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
    set_target_properties(fulldiag_evaluate PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
  endif(APPLE)
  install(TARGETS fulldiag fulldiag_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  else()
  install(TARGETS fulldiag fulldiag_evaluate RUNTIME DESTINATION bin COMPONENT applications)
  endif()
else(LAPACK_FOUND)
  message(STATUS "fulldiag will not be built since lapack library has not been found")
endif(LAPACK_FOUND)
