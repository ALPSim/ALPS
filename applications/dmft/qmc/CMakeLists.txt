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

IF(ALPS_BUILD_NGS)
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
ENDIF(ALPS_BUILD_NGS)

  set (DMFT_QMC_HYBRIDIZATION_SOURCE
       hybridization/impurity.C 
       hybridization/fourpoint.C)

IF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
  set (DMFT_QMC_HYBRIDIZATION2_SOURCE
       hybridization2/impurity.C 
       fouriertransform.C
       hybridization2/moves.C
       hybridization2/measurements.C
       hybridization2/update.C)
#  set (DMFT_QMC_HYBRIDIZATION3_SOURCE
#       hybridization3/hybsim.cpp
#       hybridization3/hyblocal.cpp
#       hybridization3/hybint.cpp
#       hybridization3/hybfun.cpp
#       hybridization3/hybmatrix.cpp
#       hybridization3/hybconfig.cpp
#       hybridization3/hybupdates.cpp
#       hybridization3/hybmeasurements.cpp)
ENDIF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)

  set (DMFT_QMC_SOLVERS "")
 
  #Interaction Expansion
  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} interaction_expansion_impl)
  add_library (interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})
  
  #Hybridization Expansion
  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} hybridization_impl)
  add_library (hybridization_impl STATIC ${DMFT_QMC_HYBRIDIZATION_SOURCE})
  
IF(ALPS_BUILD_NGS)
  #Interaction Expansion - separate solver
  add_library (interaction_expansion2_impl STATIC  ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})
ENDIF(ALPS_BUILD_NGS)
IF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
  #Hybridization Expansion - separate solver
  add_library (hybridization2_impl STATIC ${DMFT_QMC_HYBRIDIZATION2_SOURCE})
#  add_library (hybridization3_impl STATIC ${DMFT_QMC_HYBRIDIZATION3_SOURCE})
ENDIF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)

  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_definitions(${LAPACK_DEFINITIONS} )
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_executable(dmft main.C selfconsistency.C)
  add_executable(hirschfye solver_main.C hirschfyesim.C)
IF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
  add_executable(hybridization hybridization2/main.C)
#  add_executable(hybridization3 hybridization3/hybmain.cpp)
if(ALPS_LLVM_WORKAROUND)
  message(STATUS "hybridization_evaluate will not be built since llvm workaround is enabled")
else(ALPS_LLVM_WORKAROUND)
  add_executable(hybridization_evaluate
  hybridization2/hybridization_evaluate_mpi.C)
endif(ALPS_LLVM_WORKAROUND)
  add_executable(interaction interaction_expansion2/main.cpp)
ENDIF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)

  target_link_libraries(dmft dmft_qmc_impl  ${DMFT_QMC_SOLVERS})
  target_link_libraries(hirschfye dmft_qmc_impl)

IF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
  target_link_libraries(hybridization dmft_qmc_impl hybridization2_impl)
#  target_link_libraries(hybridization3 dmft_qmc_impl hybridization3_impl)
if(ALPS_LLVM_WORKAROUND)
else(ALPS_LLVM_WORKAROUND)
  target_link_libraries(hybridization_evaluate dmft_qmc_impl hybridization2_impl)
endif(ALPS_LLVM_WORKAROUND)
  target_link_libraries(interaction dmft_qmc_impl interaction_expansion2_impl)
ENDIF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)

  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
IF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
  install(TARGETS hybridization RUNTIME DESTINATION bin COMPONENT applications)
#  install(TARGETS hybridization3 RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
if(ALPS_LLVM_WORKAROUND)
else(ALPS_LLVM_WORKAROUND)
  install(TARGETS hybridization_evaluate RUNTIME DESTINATION bin COMPONENT applications)
endif(ALPS_LLVM_WORKAROUND)
ENDIF(ALPS_BUILD_NGS AND ALPS_HAVE_MPI)
#IF(ALPS_BUILD_NGS)
#  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
#ENDIF(ALPS_BUILD_NGS)

else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
