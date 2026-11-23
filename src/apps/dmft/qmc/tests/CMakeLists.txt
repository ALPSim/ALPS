# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# Regression tests for the shared DMFT QMC infrastructure.

add_executable(dmft_fouriertransform_offdiag_numeric dmft_fouriertransform_offdiag_numeric.C)
target_link_libraries(dmft_fouriertransform_offdiag_numeric PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_fouriertransform_offdiag_numeric COMPONENT dmft LABELS regression)

add_executable(dmft_hilberttransformer_itime_throw dmft_hilberttransformer_itime_throw.C)
target_link_libraries(dmft_hilberttransformer_itime_throw PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_hilberttransformer_itime_throw COMPONENT dmft LABELS regression)

add_executable(dmft_externalsolver_outfile_param dmft_externalsolver_outfile_param.C)
target_link_libraries(dmft_externalsolver_outfile_param PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_externalsolver_outfile_param COMPONENT dmft LABELS regression)

add_executable(dmft_hilberttransformer_afm_multiband_numeric
               dmft_hilberttransformer_afm_multiband_numeric.C)
target_link_libraries(dmft_hilberttransformer_afm_multiband_numeric PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_hilberttransformer_afm_multiband_numeric COMPONENT dmft LABELS regression)

add_executable(dmft_hilberttransformer_insulating dmft_hilberttransformer_insulating.C)
target_link_libraries(dmft_hilberttransformer_insulating PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_hilberttransformer_insulating COMPONENT dmft LABELS regression)

add_executable(dmft_input_files dmft_input_files.C)
target_link_libraries(dmft_input_files PRIVATE ALPS::cthyb)
alps_add_gtest(dmft_input_files COMPONENT dmft LABELS regression)

add_executable(dmft_selfconsistency_convergence dmft_selfconsistency_convergence.C
               ${CMAKE_CURRENT_SOURCE_DIR}/../selfconsistency.C)
target_link_libraries(dmft_selfconsistency_convergence PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_selfconsistency_convergence COMPONENT dmft LABELS regression)
add_executable(dmft_retarded_interaction dmft_retarded_interaction.C)
target_link_libraries(dmft_retarded_interaction PRIVATE ALPS::cthyb)
alps_add_gtest(dmft_retarded_interaction COMPONENT dmft LABELS regression)
add_executable(dmft_u_matrix dmft_u_matrix.C)
target_link_libraries(dmft_u_matrix PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_u_matrix COMPONENT dmft LABELS regression)

add_executable(dmft_hilberttransformer_uniform_field dmft_hilberttransformer_uniform_field.C)
target_link_libraries(dmft_hilberttransformer_uniform_field PRIVATE dmft_qmc_impl)
alps_add_gtest(dmft_hilberttransformer_uniform_field COMPONENT dmft LABELS regression)

# dmft stops when SYMMETRIZATION is on together with a field.
foreach(field H H_INIT)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dmft_symmetrized_${field}/parm"
    "SOLVER = \"Hirsch-Fye\"\nOMEGA_LOOP = 1\nSYMMETRIZATION = true\nMAX_IT = 1\nCONVERGED = 0.001\n"
    "BETA = 4\nU = 2\nt = 0.5\nMU = 0\nN = 16\nNMATSUBARA = 40\nSWEEPS = 100\nTHERMALIZATION = 10\n"
    "${field} = 0.05\n")
  add_test(NAME dmft_rejects_symmetrized_${field}
           COMMAND dmft "${CMAKE_CURRENT_BINARY_DIR}/dmft_symmetrized_${field}/parm")
  set_tests_properties(dmft_rejects_symmetrized_${field} PROPERTIES LABELS "dmft;regression" TIMEOUT 60
    PASS_REGULAR_EXPRESSION "SYMMETRIZATION = true requires H = 0 and H_INIT = 0")
endforeach()

add_executable(dmft_interaction_expansion_choice dmft_interaction_expansion_choice.C)
alps_add_gtest(dmft_interaction_expansion_choice COMPONENT dmft LABELS regression)

add_executable(dmft_ctint_sz_results dmft_ctint_sz_results.C)
target_link_libraries(dmft_ctint_sz_results PRIVATE ALPS::ctint)
alps_add_gtest(dmft_ctint_sz_results COMPONENT dmft LABELS regression)

add_executable(dmft_ctint_max_order dmft_ctint_max_order.C)
target_link_libraries(dmft_ctint_max_order PRIVATE ALPS::ctint)
if(NOT ALPS_ENABLE_MPI)
  alps_add_gtest(dmft_ctint_max_order COMPONENT dmft LABELS regression)
else()
  alps_add_mpi_gtest(dmft_ctint_max_order COMPONENT dmft RANKS 2 LABELS regression)
endif()

# The serial solver stops after SWEEPS; under MPI it stops at the adapter's
# scheduled checks, where a rank can still be thermalizing.
if(NOT ALPS_ENABLE_MPI)
  add_executable(dmft_ctint_atomic dmft_ctint_atomic.C)
  target_link_libraries(dmft_ctint_atomic PRIVATE ALPS::ctint)
  alps_add_gtest(dmft_ctint_atomic COMPONENT dmft LABELS regression)

  add_executable(dmft_ctint_thermalization dmft_ctint_thermalization.C)
  target_link_libraries(dmft_ctint_thermalization PRIVATE ALPS::ctint)
  alps_add_gtest(dmft_ctint_thermalization COMPONENT dmft LABELS regression)
else()
  add_executable(dmft_ctint_thermalization_mpi dmft_ctint_thermalization_mpi.C)
  target_link_libraries(dmft_ctint_thermalization_mpi PRIVATE ALPS::ctint)
  alps_add_mpi_gtest(dmft_ctint_thermalization_mpi COMPONENT dmft RANKS 2 LABELS regression)
endif()

if(TARGET ALPS::cthyb_segment)
  add_executable(dmft_cthyb_segment_exact dmft_cthyb_segment_exact.C)
  target_link_libraries(dmft_cthyb_segment_exact PRIVATE ALPS::cthyb_segment)
  if(NOT ALPS_ENABLE_MPI)
    alps_add_gtest(dmft_cthyb_segment_exact COMPONENT dmft LABELS regression TIMEOUT 300)
  else()
    alps_add_mpi_gtest(dmft_cthyb_segment_exact COMPONENT dmft RANKS 2 LABELS regression TIMEOUT 300)
  endif()
endif()
