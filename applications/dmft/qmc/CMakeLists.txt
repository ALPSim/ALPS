#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)

  set (DMFT_QMC_SOURCE
       externalsolver.C 
       fouriertransform.C 
       hilberttransformer.C 
       hirschfyesim.C 
       alps_solver.C 
       bandstructure.C 
       auxiliaryfunctions.C)
  add_library(dmft_qmc_impl STATIC ${DMFT_QMC_SOURCE})

  set (DMFT_QMC_INTERACTION_EXPANSION_SOURCE 
       interaction_expansion/auxiliary.cpp
       interaction_expansion/observables.cpp 
       interaction_expansion/fastupdate.cpp
       interaction_expansion/selfenergy.cpp 
       interaction_expansion/solver.cpp 
       interaction_expansion/io.cpp 
       interaction_expansion/splines.cpp 
       interaction_expansion/interaction_expansion.cpp 
       interaction_expansion/measurements.cpp 
      interaction_expansion/model.cpp)

IF(ALPS_HAVE_MPI)
  set (DMFT_QMC_INTERACTION_EXPANSION2_SOURCE
       fouriertransform.C
       interaction_expansion2/auxiliary.cpp
       interaction_expansion2/observables.cpp 
       interaction_expansion2/fastupdate.cpp
       interaction_expansion2/selfenergy.cpp 
       interaction_expansion2/solver.cpp 
       interaction_expansion2/io.cpp 
       interaction_expansion2/splines.cpp 
       interaction_expansion2/interaction_expansion.cpp 
       interaction_expansion2/measurements.cpp 
       interaction_expansion2/model.cpp)
ENDIF(ALPS_HAVE_MPI)

  set (DMFT_QMC_HYBRIDIZATION_SOURCE
       hybridization/hybsim.cpp
       hybridization/hyblocal.cpp
       hybridization/hybint.cpp
       hybridization/hybfun.cpp
       hybridization/hybretintfun.cpp
       hybridization/hybmatrix.cpp
       hybridization/hybmatrix_ft.cpp
       hybridization/hybconfig.cpp
       hybridization/hybupdates.cpp
       hybridization/hybevaluate.cpp
       hybridization/hybmeasurements.cpp)

  set (DMFT_QMC_SOLVERS "")
 
  #Interaction Expansion
  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} interaction_expansion_impl)
  add_library (interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})

  #Interaction Expansion - separate solver
IF(ALPS_HAVE_MPI)
  add_library (interaction_expansion2_impl STATIC  ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})
  IF(ALPS_BUILD_PYTHON)
    add_library (ctint.so SHARED interaction_expansion2/main.cpp ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})
    SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
    set_target_properties(ctint.so PROPERTIES PREFIX "")
    set_target_properties(ctint.so PROPERTIES SUFFIX "")
    set_target_properties(ctint.so PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE")
    target_link_libraries(ctint.so alps ${ALPS_LIBRARIES} ${ALPS_BLAS_LIBRARY} ${BLAS_LIBRARY} ${ALPS_LAPACK_LIBRARY} ${LAPACK_LIBRARY})
  ENDIF(ALPS_BUILD_PYTHON)
ENDIF(ALPS_HAVE_MPI)

IF(CMAKE_VERSION VERSION_GREATER "2.8.8")
      #Hybridization Expansion - separate solver
      add_library (hybridization_impl_obj OBJECT ${DMFT_QMC_HYBRIDIZATION_SOURCE})
      set_property(TARGET hybridization_impl_obj PROPERTY POSITION_INDEPENDENT_CODE 1)
      add_library (hybridization_impl STATIC $<TARGET_OBJECTS:hybridization_impl_obj>)
      #Hybridization Expansion - python module
      IF(ALPS_BUILD_PYTHON)
        add_library (cthyb.so SHARED hybridization/hybmain.cpp $<TARGET_OBJECTS:hybridization_impl_obj>)
        SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
        set_target_properties(cthyb.so PROPERTIES PREFIX "")
        set_target_properties(cthyb.so PROPERTIES SUFFIX "")
        set_target_properties(cthyb.so PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE")
        target_link_libraries(cthyb.so alps ${ALPS_LIBRARIES} ${ALPS_BLAS_LIBRARY} ${BLAS_LIBRARY} ${ALPS_LAPACK_LIBRARY} ${LAPACK_LIBRARY})
      ENDIF(ALPS_BUILD_PYTHON)
ELSE(CMAKE_VERSION VERSION_GREATER "2.8.8")
    IF(ALPS_HAVE_MPI)
      #Hybridization Expansion - separate solver
      add_library (hybridization_impl STATIC ${DMFT_QMC_HYBRIDIZATION_SOURCE})
      #Hybridization Expansion - python module
      IF(ALPS_BUILD_PYTHON)
        add_library (cthyb.so SHARED hybridization/hybmain.cpp ${DMFT_QMC_HYBRIDIZATION_SOURCE})
        SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
        set_target_properties(cthyb.so PROPERTIES PREFIX "")
        set_target_properties(cthyb.so PROPERTIES SUFFIX "")
        set_target_properties(cthyb.so PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE")
        target_link_libraries(cthyb.so alps ${ALPS_LIBRARIES} ${ALPS_BLAS_LIBRARY} ${BLAS_LIBRARY} ${ALPS_LAPACK_LIBRARY} ${LAPACK_LIBRARY})
      ENDIF(ALPS_BUILD_PYTHON)
    ENDIF(ALPS_HAVE_MPI)
ENDIF(CMAKE_VERSION VERSION_GREATER "2.8.8")

  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_definitions(${LAPACK_DEFINITIONS} )
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_executable(dmft main.C selfconsistency.C)
  add_executable(hirschfye solver_main.C hirschfyesim.C)

IF(ALPS_HAVE_MPI)
  add_executable(hybridization hybridization/hybmain.cpp)
  add_executable(interaction interaction_expansion2/main.cpp)
ENDIF(ALPS_HAVE_MPI)

  target_link_libraries(dmft dmft_qmc_impl  ${DMFT_QMC_SOLVERS})
  target_link_libraries(hirschfye dmft_qmc_impl)

IF(ALPS_HAVE_MPI)
  target_link_libraries(hybridization dmft_qmc_impl hybridization_impl)
  target_link_libraries(interaction dmft_qmc_impl interaction_expansion2_impl)
ENDIF(ALPS_HAVE_MPI)

  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
IF(ALPS_HAVE_MPI)
  install(TARGETS hybridization RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
  IF(ALPS_BUILD_PYTHON)
    install(TARGETS cthyb.so LIBRARY DESTINATION lib/pyalps COMPONENT applications)
    install(TARGETS ctint.so LIBRARY DESTINATION lib/pyalps COMPONENT applications)
  ENDIF(ALPS_BUILD_PYTHON)
ENDIF(ALPS_HAVE_MPI)
  INSTALL(FILES hybridization/Documentation/hybdoc.pdf DESTINATION doc)
else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
