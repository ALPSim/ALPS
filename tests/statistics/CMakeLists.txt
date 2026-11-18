# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# Statistics against exact values, on the same input for every statistics
# library in the tree: ALEA, alps::accumulators and ALEA v2.
if(NOT TARGET ALPS::aleav2)
  return()
endif()

foreach(name exact_moments same_stream ar1_known_answers coverage nonlinear vector_covariance
             merge signed_ratio hdf5_roundtrip frozen_files observables observables_files)
  add_executable(statistics_${name} ${name}.cpp)
  target_link_libraries(statistics_${name} PRIVATE ALPS::alps ALPS::aleav2)
  alps_add_gtest(statistics_${name} COMPONENT statistics TIMEOUT 120)
endforeach()

# Wrote fixtures/legacy_observables.{h5,expected}; built to keep it compiling,
# never run by a test.
add_executable(statistics_write_legacy_observables fixtures/write_legacy_observables.cpp)
target_link_libraries(statistics_write_legacy_observables PRIVATE ALPS::alps)
target_compile_features(statistics_write_legacy_observables PRIVATE cxx_std_17)

if(ALPS_ENABLE_MPI)
  add_executable(statistics_mpi_reduce mpi_reduce.cpp)
  target_link_libraries(statistics_mpi_reduce PRIVATE ALPS::alps ALPS::aleav2)
  alps_add_mpi_gtest(statistics_mpi_reduce COMPONENT statistics RANKS 2)
  alps_add_mpi_gtest(statistics_mpi_reduce COMPONENT statistics RANKS 3)
endif()
