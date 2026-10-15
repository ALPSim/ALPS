#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
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

foreach(name default ordering stream assign not_found external checkpoint)
  add_executable(param_${name} ${name}.cpp)

  target_link_libraries(param_${name} PRIVATE ALPS::params)
  alps_add_gtest(param_${name} COMPONENT params)
endforeach()

# ALPSCore's params tests, copied unchanged.
foreach(name serialize_variant dictionary params dots_in_keys params_default_vectors
             params_string_vectors params_formats params_cmdline params_defaulted
             params_has_unused params_origins convenience_params dictval_eq dictionary_eq
             params_eq dictionary_hdf5 params_hdf5 params_hdf5_with_ini params_vec_hdf5)
  if(name MATCHES "^params")
    set(target ${name})
  else()
    set(target params_${name})
  endif()
  add_executable(${target} ${name}.cpp)
  target_link_libraries(${target} PRIVATE ALPS::params)
  # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
  target_compile_options(${target} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  # ALPSCore's alps::testing::unique_file, from tests/support, which the
  # tests' params_test_support.hpp uses.
  target_sources(${target} PRIVATE "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")
  alps_add_gtest(${target} COMPONENT params)
endforeach()

# ALPS's own tests of its additions to ALPSCore's params.
foreach(name params_supply params_plain_group)
  add_executable(${name} ${name}.cpp "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")
  target_link_libraries(${name} PRIVATE ALPS::params)
  alps_add_gtest(${name} COMPONENT params)
endforeach()

if(ALPS_ENABLE_MPI)
  # They define main() themselves, as ALPSCore's MPI tests do.
  foreach(name serialize_variant_mpi dictionary_bcast_mpi params_bcast_mpi)
    if(name MATCHES "^params")
      set(target ${name})
    else()
      set(target params_${name})
    endif()
    add_executable(${target} ${name}.cpp "${PROJECT_SOURCE_DIR}/tests/support/src/unique_file.cpp")
    target_link_libraries(${target} PRIVATE ALPS::params)
    alps_add_mpi_gtest(${target} COMPONENT params RANKS 2 OWN_MAIN)
    alps_add_mpi_gtest(${target} COMPONENT params RANKS 3 OWN_MAIN)
  endforeach()
endif()
