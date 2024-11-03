#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

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
  set_target_properties(spinmc  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  set_target_properties(spinmc_evaluate  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  install(TARGETS spinmc spinmc_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  else()
  install(TARGETS spinmc spinmc_evaluate RUNTIME DESTINATION bin COMPONENT applications)
  endif()
else(LAPACK_FOUND)
  message(STATUS "spins will not be built since lapack library is not found")
endif(LAPACK_FOUND)
