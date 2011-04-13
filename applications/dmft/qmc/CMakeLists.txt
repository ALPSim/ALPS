#  cOPYRIGHT mATTHIAS tROYER 2009 - 2010.
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
       2dsimpson.C 
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
  set (DMFT_QMC_INTERACTION_EXPANSION2_SOURCE 
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
  set (DMFT_QMC_HYBRIDIZATION_SOURCE
       hybridization/impurity.C 
       hybridization/fourpoint.C)
  set (DMFT_QMC_HYBRIDIZATION2_SOURCE
       hybridization2/impurity.C 
       fouriertransform.C
       hybridization2/moves.C
       hybridization2/update.C)

  set (DMFT_QMC_SOLVERS "")
 
  #Interaction Expansion
  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} interaction_expansion_impl)
  add_library (interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})
  
  #Hybridization Expansion
  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} hybridization_impl)
  add_library (hybridization_impl STATIC ${DMFT_QMC_HYBRIDIZATION_SOURCE})
  
  #Interaction Expansion - separate solver
  add_library (interaction_expansion2_impl STATIC  ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})

  #Hybridization Expansion - separate solver
  add_library (hybridization2_impl STATIC ${DMFT_QMC_HYBRIDIZATION2_SOURCE})

  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_definitions(${LAPACK_DEFINITIONS} )
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_executable(dmft main.C selfconsistency.C)
  add_executable(hirschfye solver_main.C hirschfyesim.C)
  add_executable(hybridization hybridization2/main.C)
  add_executable(interaction interaction_expansion2/main.cpp)

  target_link_libraries(dmft dmft_qmc_impl  ${DMFT_QMC_SOLVERS})
  target_link_libraries(hirschfye dmft_qmc_impl)
  target_link_libraries(hybridization dmft_qmc_impl hybridization2_impl)
  target_link_libraries(interaction dmft_qmc_impl interaction_expansion2_impl)

  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hybridization RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
