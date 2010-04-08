#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)

  add_library(dmft_qmc_impl STATIC externalsolver.C fouriertransform.C hilberttransformer.C hirschfyesim.C alps_solver.C 2dsimpson.C auxiliaryfunctions.C)
  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

  add_definitions(${LAPACK_DEFINITIONS} -DBUILD_DMFT_QMC_HYBRIDIZATION)
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_executable(dmft main.C selfconsistency.C hybridization/impurity.C hybridization/fourpoint.C)
add_executable(hirschfye solver_main.C hirschfyesim.C)

  target_link_libraries(dmft dmft_qmc_impl)
  target_link_libraries(hirschfye dmft_qmc_impl)
  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
