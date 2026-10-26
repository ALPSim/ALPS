# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's accumulators tests, copied with publication pointers removed, and ALPS's own
# ar1_calibration and linear_transform_error. The legacy_alea/ and legacy_accumulator/ directories hold
# older NGS tests and are not built.
foreach(name testharness assign_acc get_results add_accumulators count divide_accumulators
             unary_ops mean merge mult_by_constant weighted nonexistent acc_count_test
             mean_err_count save_load save_load2 vec_const_binop_simple binop_with_constant
             binop_with_scalar binop_mixed binop_mixed_faildemo single_accumulator
             autocorrelation concurrent_access print scalar_result_type negative_error
             ar1_calibration linear_transform_error)
  add_executable(accumulators_${name} ${name}.cpp)
  target_link_libraries(accumulators_${name} PRIVATE ALPS::accumulators)
  # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
  target_compile_options(accumulators_${name} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  if(name STREQUAL "save_load2")
    # Its typed cases write the same file names in the working directory.
    # ALPSCore ran them in one process; CTest runs each case as its own test.
    alps_add_gtest(accumulators_${name} COMPONENT accumulators
      PROPERTIES RESOURCE_LOCK accumulators_save_load2)
  else()
    alps_add_gtest(accumulators_${name} COMPONENT accumulators)
  endif()
endforeach()
# ALPSCore's alps::testing::unique_file, from tests/support.
target_sources(accumulators_save_load PRIVATE
  "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")

if(ALPS_ENABLE_MPI)
  # They define main() themselves, as ALPSCore's MPI tests do.
  foreach(name mpi_merge mpi_merge_uneven repeated_merge zero_vector_mpi)
    add_executable(accumulators_${name} ${name}.cpp)
    target_link_libraries(accumulators_${name} PRIVATE ALPS::accumulators)
    target_compile_options(accumulators_${name} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
    alps_add_mpi_gtest(accumulators_${name} COMPONENT accumulators RANKS 2 OWN_MAIN)
    alps_add_mpi_gtest(accumulators_${name} COMPONENT accumulators RANKS 3 OWN_MAIN)
  endforeach()
endif()
