# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's ALEA tests, plus local regressions.
foreach(name compare_test fisher_f empty twogauss model result transform
             stream_serializer)
  add_executable(alea_${name} ${name}.cpp)
  target_link_libraries(alea_${name} PRIVATE ALPS::alea)
  # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
  target_compile_options(alea_${name} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  if(name STREQUAL "twogauss")
    # Every twogauss case writes twogauss.hdf5 in the working directory. ALPSCore
    # ran the cases in one process; CTest runs each case as its own test.
    alps_add_gtest(alea_${name} COMPONENT alea PROPERTIES RESOURCE_LOCK alea_twogauss_hdf5)
  else()
    alps_add_gtest(alea_${name} COMPONENT alea)
  endif()
endforeach()

if(ALPS_ENABLE_MPI)
  # It defines main() itself, as ALPSCore's MPI tests do.
  add_executable(alea_mpi_twogauss mpi_twogauss.cpp)
  target_link_libraries(alea_mpi_twogauss PRIVATE ALPS::alea)
  target_compile_options(alea_mpi_twogauss PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  alps_add_mpi_gtest(alea_mpi_twogauss COMPONENT alea RANKS 2 OWN_MAIN)
  alps_add_mpi_gtest(alea_mpi_twogauss COMPONENT alea RANKS 3 OWN_MAIN)
endif()
