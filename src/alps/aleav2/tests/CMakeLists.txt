# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's ALEA tests, copied unchanged, plus coexistence with ALPS's ALEA.
foreach(name compare_test fisher_f empty twogauss galois model result transform
             stream_serializer)
  add_executable(aleav2_${name} ${name}.cpp)
  target_link_libraries(aleav2_${name} PRIVATE ALPS::aleav2)
  # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
  target_compile_options(aleav2_${name} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  if(name STREQUAL "twogauss")
    # Every twogauss case writes twogauss.hdf5 in the working directory. ALPSCore
    # ran the cases in one process; CTest runs each case as its own test.
    alps_add_gtest(aleav2_${name} COMPONENT aleav2 PROPERTIES RESOURCE_LOCK aleav2_twogauss_hdf5)
  else()
    alps_add_gtest(aleav2_${name} COMPONENT aleav2)
  endif()
endforeach()

add_executable(aleav2_coexistence coexistence.cpp)
target_link_libraries(aleav2_coexistence PRIVATE ALPS::aleav2 ALPS::alps)
alps_add_gtest(aleav2_coexistence COMPONENT aleav2)

if(ALPS_ENABLE_MPI)
  # It defines main() itself, as ALPSCore's MPI tests do.
  add_executable(aleav2_mpi_twogauss mpi_twogauss.cpp)
  target_link_libraries(aleav2_mpi_twogauss PRIVATE ALPS::aleav2)
  target_compile_options(aleav2_mpi_twogauss PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  alps_add_mpi_gtest(aleav2_mpi_twogauss COMPONENT aleav2 RANKS 2 OWN_MAIN)
  alps_add_mpi_gtest(aleav2_mpi_twogauss COMPONENT aleav2 RANKS 3 OWN_MAIN)
endif()
