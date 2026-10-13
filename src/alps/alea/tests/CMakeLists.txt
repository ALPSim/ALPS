# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPSCore's ALEA tests, copied unchanged, plus coexistence with legacy ALEA.
# ALPSCore's mpi_twogauss is not ported: ALPS does not build alea/mpi.hpp.
foreach(name compare_test fisher_f empty twogauss galois model result transform
             stream_serializer)
  add_executable(alea_${name} ${name}.cpp)
  target_link_libraries(alea_${name} PRIVATE ALPS::alea)
  # alps/testing/near.hpp, ALPSCore's test-only comparison helpers.
  target_include_directories(alea_${name} PRIVATE support)
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

add_executable(alea_legacy_coexistence legacy_coexistence.cpp)
target_link_libraries(alea_legacy_coexistence PRIVATE ALPS::alea ALPS::alps)
alps_add_gtest(alea_legacy_coexistence COMPONENT alea)
