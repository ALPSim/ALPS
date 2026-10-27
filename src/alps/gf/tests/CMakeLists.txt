# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's gf tests, copied with publication pointers removed. support/ holds
# ALPSCore's unique_file helper, also unchanged; the helpers other modules share
# are in the top-level tests/support.
alps_require_googletest()

# ALPSCore compiled gf_test.cpp and its helpers into every test; so does this.
add_library(alps_gf_test_support OBJECT gf_test.cpp support/src/unique_file.cpp)
target_include_directories(alps_gf_test_support PUBLIC support "${PROJECT_SOURCE_DIR}/tests/support")
target_link_libraries(alps_gf_test_support PUBLIC ALPS::gf GTest::gtest)
target_compile_definitions(alps_gf_test_support PUBLIC ${ALPS_GF_ASSERT_DEFINITION})
# gf1_print_test keeps ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
target_compile_options(alps_gf_test_support PUBLIC
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)

# gf_check_eigen has no test cases: it only checks that Eigen's headers compile.
# ALPSCore ran it as an empty test; CTest would find nothing to run.
add_library(gf_check_eigen OBJECT gf_check_eigen.cpp)
target_link_libraries(gf_check_eigen PRIVATE alps_gf_test_support)

foreach(name gf_new_test gf_new_tail_test index_test mesh_test one_index_gf_test
             gf1_print_test two_index_gf_test three_index_gf_test four_index_gf_test
             four_index_gf_pn_test seven_index_gf_test itime_gf_test fourier_test
             grid_test piecewise_polynomial_test)
  # Prefix the target names, which are global, unless they already start with gf.
  if(name MATCHES "^gf")
    set(target ${name})
  else()
    set(target gf_${name})
  endif()
  add_executable(${target} ${name}.cpp)
  target_link_libraries(${target} PRIVATE alps_gf_test_support)
  # CTest runs each case as its own test. Cases that write a fixed file name each
  # use their own, so they can run in parallel.
  alps_add_gtest(${target} COMPONENT gf)
endforeach()

if(ALPS_ENABLE_MPI)
  # The MPI tests define main() and count broadcasts by intercepting MPI_Init and
  # MPI_Bcast in mpi_guard.cpp, so its object file goes into each test.
  add_library(alps_gf_mpi_test_support OBJECT mpi_guard.cpp)
  target_link_libraries(alps_gf_mpi_test_support PUBLIC alps_gf_test_support MPI::MPI_CXX)
  foreach(name one_index_gf_test_mpi four_index_gf_test_mpi four_index_gf_test_mismatched_mpi
               four_index_gf_test_mismatched-tail_mpi multiarray_bcast_mpi mesh_test_mpi
               gf_new_test_mpi gf_new_tail_test_mpi)
    if(name MATCHES "^gf")
      set(target ${name})
    else()
      set(target gf_${name})
    endif()
    add_executable(${target} ${name}.cpp)
    # Object files reach only the targets that link their library directly.
    target_link_libraries(${target} PRIVATE alps_gf_test_support alps_gf_mpi_test_support)
    alps_add_mpi_gtest(${target} COMPONENT gf RANKS 2 OWN_MAIN)
    alps_add_mpi_gtest(${target} COMPONENT gf RANKS 3 OWN_MAIN)
  endforeach()

endif()
