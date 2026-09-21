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

set(DMFT_QMC_SOURCE
    externalsolver.C
    fouriertransform.C
    hilberttransformer.C
    hirschfyesim.C
    alps_solver.C
    bandstructure.C
    auxiliaryfunctions.C)
add_library(dmft_qmc_impl STATIC ${DMFT_QMC_SOURCE})

set(DMFT_QMC_INTERACTION_EXPANSION_SOURCE
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

add_library(alps_cthyb STATIC
    hybridization/run.cpp
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
add_library(alps_ctint STATIC
    interaction_expansion2/run.cpp
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
foreach(solver IN ITEMS cthyb ctint)
  add_library(ALPS::${solver} ALIAS alps_${solver})
  set_target_properties(alps_${solver} PROPERTIES EXPORT_NAME ${solver}
    CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
  target_link_libraries(alps_${solver} PUBLIC ALPS::alps)
  target_include_directories(alps_${solver} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
endforeach()
install(TARGETS alps_cthyb alps_ctint EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

#Interaction Expansion
add_library(interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})
target_link_libraries(interaction_expansion_impl PUBLIC ALPS::alps)

target_link_libraries(dmft_qmc_impl PUBLIC ALPS::alps)

target_include_directories(dmft_qmc_impl PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(interaction_expansion_impl PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(dmft main.C selfconsistency.C)
add_executable(hirschfye solver_main.C)

add_executable(hybridization hybridization/hybmain.cpp)
add_executable(interaction interaction_expansion2/main.cpp)

target_link_libraries(dmft PRIVATE dmft_qmc_impl interaction_expansion_impl)
target_link_libraries(hirschfye PRIVATE dmft_qmc_impl)

target_link_libraries(hybridization PRIVATE ALPS::cthyb)
target_link_libraries(interaction PRIVATE ALPS::ctint)

if(BUILD_TESTING)
  # Regression tests for the shared DMFT QMC infrastructure.
  add_executable(dmft_fouriertransform_offdiag_numeric dmft_fouriertransform_offdiag_numeric.C)
  target_link_libraries(dmft_fouriertransform_offdiag_numeric PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_fouriertransform_offdiag_numeric)
  set_property(TEST dmft_fouriertransform_offdiag_numeric PROPERTY LABELS dmft)

  add_executable(dmft_hilberttransformer_itime_throw dmft_hilberttransformer_itime_throw.C)
  target_link_libraries(dmft_hilberttransformer_itime_throw PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_hilberttransformer_itime_throw)
  set_property(TEST dmft_hilberttransformer_itime_throw PROPERTY LABELS dmft)

  add_executable(dmft_externalsolver_outfile_param dmft_externalsolver_outfile_param.C)
  target_link_libraries(dmft_externalsolver_outfile_param PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_externalsolver_outfile_param)
  set_property(TEST dmft_externalsolver_outfile_param PROPERTY LABELS dmft)

  add_executable(dmft_hilberttransformer_afm_multiband_numeric
                 dmft_hilberttransformer_afm_multiband_numeric.C)
  target_link_libraries(dmft_hilberttransformer_afm_multiband_numeric PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_hilberttransformer_afm_multiband_numeric)
  set_property(TEST dmft_hilberttransformer_afm_multiband_numeric PROPERTY LABELS dmft)

  add_executable(dmft_interaction_expansion_choice dmft_interaction_expansion_choice.C)
  alps_add_test(dmft_interaction_expansion_choice)
  set_property(TEST dmft_interaction_expansion_choice PROPERTY LABELS dmft)
endif()
install(FILES hybridization/Documentation/hybdoc.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})

alps_install_applications(dmft hirschfye hybridization interaction)
