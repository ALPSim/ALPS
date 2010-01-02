if(LAPACK_FOUND)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(spinmc spinmc.C abstract_fitter.C base_incr_fitter.C dummy_fitter.C
		 est_grad_fitter.C fitting_scheduler.C spinmc_factory.C)
  target_link_libraries(spinmc alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_executable(spinmc_evaluate spinmc_evaluate.C abstract_fitter.C base_incr_fitter.C
		 dummy_fitter.C est_grad_fitter.C fitting_scheduler.C
		 spinmc_factory.C)
  target_link_libraries(spinmc_evaluate alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  install(TARGETS spinmc spinmc_evaluate RUNTIME DESTINATION bin)
else(LAPACK_FOUND)
  message(STATUS "spins will not be built since lapack library is not found")
endif(LAPACK_FOUND)
