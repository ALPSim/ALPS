if(LAPACK_FOUND)
  set (WERNER TRUE)
  set (RUBTSOV TRUE)

  add_subdirectory(solvers)
  add_definitions(${LAPACK_DEFINITIONS})
  include_directories(${PROJECT_SOURCE_DIR}/applications/dca)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS} -L/opt/local/lib -lgsl")

  set (WERNER_SOURCE solvers/Werner/impurity.C solvers/Werner/fourpoint.C)
  set (RUBTSOV_SOURCE 
      solvers/Rubtsov/rubtsov_auxiliary.cpp
      solvers/Rubtsov/rubtsov_observables.cpp 
      solvers/Rubtsov/rubtsov_fastupdate.cpp
      solvers/Rubtsov/rubtsov_fourpoint.cpp 
      solvers/Rubtsov/rubtsov_selfenergy.cpp 
      solvers/Rubtsov/rubtsov_solver.cpp 
      solvers/Rubtsov/rubtsov_io.cpp 
      solvers/Rubtsov/rubtsov_splines.cpp 
      solvers/Rubtsov/rubtsov.cpp 
      solvers/Rubtsov/rubtsov_measurements.cpp 
      solvers/Rubtsov/rubtsov_wang_landau.cpp 
      solvers/Rubtsov/rubtsov_model.cpp 
      solvers/Rubtsov/rubtsov_worm.cpp)
  set (DCA_SOURCE 2dsimpson.C selfconsistency.C clusterhelper.C
       	clustertransformer.C externalsolver.C hirschfyesim.C
  	alps_solver.C fouriertransform.C
	auxiliaryfunctions.C hilberttransformer.C)
  add_library(dca_impl STATIC ${DCA_SOURCE})
  set (DCA_SOLVERS "")
  if(RUBTSOV)
    set (DCA_SOLVERS ${DCA_SOLVERS} rubtsov_impl)
    add_library (rubtsov_impl STATIC ${RUBTSOV_SOURCE})
  endif(RUBTSOV)
  if(WERNER)
    set (DCA_SOLVERS ${DCA_SOLVERS} werner_impl)
    add_library (werner_impl STATIC ${WERNER_SOURCE})
  endif(WERNER)
  add_executable(dca main.C)
  add_executable(interpolate interpolate.C interpolator.C)
  add_executable(generate generate.C generator2D.C generator3D.C)
  target_link_libraries(dca dca_impl ${DCA_SOLVERS} alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  target_link_libraries(interpolate dca_impl alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  target_link_libraries(generate dca_impl alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  install(TARGETS dca RUNTIME DESTINATION bin)
  install(TARGETS interpolate RUNTIME DESTINATION bin)
  install(TARGETS generate RUNTIME DESTINATION bin)
  install(FILES clusters.xml DESTINATION lib/xml)
else(LAPACK_FOUND)
  message(STATUS "dca will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
