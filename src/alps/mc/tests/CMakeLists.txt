# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's mc tests, copied unchanged.
foreach(name sum_single timer_in_sim timer check_schedule)
  add_executable(mc_${name} ${name}.cpp)
  target_link_libraries(mc_${name} PRIVATE ALPS::mc)
  # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
  target_compile_options(mc_${name} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  alps_add_gtest(mc_${name} COMPONENT mc)
endforeach()
# ALPSCore's alps::testing::unique_file, from tests/support.
target_sources(mc_sum_single PRIVATE "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")

if(ALPS_ENABLE_MPI)
  # They define main() themselves, as ALPSCore's MPI tests do. reduce and
  # signed_obs have no GoogleTest cases; their exit status is the check.
  foreach(name sum_mpi reduce signed_obs custom_scheduler reduce_unavailable_results)
    add_executable(mc_${name} ${name}.cpp "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")
    target_link_libraries(mc_${name} PRIVATE ALPS::mc)
    alps_add_mpi_gtest(mc_${name} COMPONENT mc RANKS 2 OWN_MAIN)
    alps_add_mpi_gtest(mc_${name} COMPONENT mc RANKS 3 OWN_MAIN)
  endforeach()
endif()
