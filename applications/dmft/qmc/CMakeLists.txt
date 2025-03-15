#  Copyright Matthias Troyer 2009 - 2013.
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
  target_link_libraries(interaction_expansion_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  add_definitions(${LAPACK_DEFINITIONS} )
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_executable(dmft main.C selfconsistency.C)
  add_executable(hirschfye solver_main.C hirschfyesim.C)

  add_executable(hybridization hybridization/hybmain.cpp ${DMFT_QMC_HYBRIDIZATION_SOURCE})
  add_executable(interaction interaction_expansion2/main.cpp ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})

  target_link_libraries(dmft dmft_qmc_impl  ${DMFT_QMC_SOLVERS})
  target_link_libraries(hirschfye dmft_qmc_impl)

  target_link_libraries(hybridization dmft_qmc_impl)
  target_link_libraries(interaction dmft_qmc_impl)
if(ALPS_PYTHON_WHEEL)
  set_target_properties(dmft  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_BINARY_DIR}/src/boost" )
  set_target_properties(hirschfye  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  set_target_properties(hybridization  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  set_target_properties(interaction  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  install(TARGETS dmft RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  install(TARGETS hybridization RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  install(TARGETS interaction RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS hybridization RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
  install(FILES hybridization/Documentation/hybdoc.pdf DESTINATION doc)
endif()
else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
