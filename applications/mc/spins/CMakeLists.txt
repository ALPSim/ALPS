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
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_library(spinmc_impl STATIC abstract_fitter.C base_incr_fitter.C
		 dummy_fitter.C est_grad_fitter.C fitting_scheduler.C spinmc_factory.C)
  target_link_libraries(spinmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_executable(spinmc spinmc.C)
  add_executable(spinmc_evaluate spinmc_evaluate.C)
  target_link_libraries(spinmc spinmc_impl)
  target_link_libraries(spinmc_evaluate spinmc_impl)
  if(ALPS_PYTHON_WHEEL)
  set_target_properties(spinmc  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
  set_target_properties(spinmc_evaluate  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
  install(TARGETS spinmc spinmc_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  else()
  install(TARGETS spinmc spinmc_evaluate RUNTIME DESTINATION bin COMPONENT applications)
  endif()
else(LAPACK_FOUND)
  message(STATUS "spins will not be built since lapack library is not found")
endif(LAPACK_FOUND)
