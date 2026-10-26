#  Copyright Synge Todo 2012.
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

foreach(name vmusage)
  add_executable(${name} ${name}.cpp)

  target_link_libraries(${name} PRIVATE ALPS::utilities)
  alps_add_gtest(${name} COMPONENT utilities)
endforeach()

add_executable(ngs_hash ngs_hash.cpp)
target_link_libraries(ngs_hash PRIVATE ALPS::utilities)
alps_add_gtest(ngs_hash COMPONENT utilities)

# This consumer has no aggregate or numerical include roots.
add_executable(foundation_headers foundation_headers.cpp)
target_link_libraries(foundation_headers PRIVATE ALPS::utilities)
alps_add_gtest(foundation_headers COMPONENT utilities LABELS architecture)

# Link only the owning foundation: no aggregate library may hide missing symbols.
add_executable(citations citations.cpp)
target_link_libraries(citations PRIVATE ALPS::utilities)
alps_add_gtest(citations COMPONENT utilities LABELS citations)

if(ALPS_ENABLE_MPI)
  # ALPSCore's tests of alps/utilities/mpi.hpp, copied with publication pointers removed with their
  # helpers. They define main() themselves; the environment tests start and stop
  # MPI inside their single case. mpi_utils_bcast_optional is not ported: it
  # tests ALPSCore's alps/utilities/mpi_optional.hpp, which ALPS does not have.
  foreach(name mpi_utils_env mpi_utils_comm mpi_utils_bcast mpi_utils_reduce
               mpi_utils_env_exceptions1 mpi_utils_env_exceptions2)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE ALPS::utilities)
    if(name MATCHES "exceptions")
      # mpi_test_support.cpp intercepts MPI_Abort to count the calls.
      target_sources(${name} PRIVATE mpi_test_support.cpp)
    endif()
    if(name STREQUAL "mpi_utils_env_exceptions1")
      # The only one without a main(), as in ALPSCore.
      target_link_libraries(${name} PRIVATE GTest::gtest_main)
    endif()
    # The copies keep ALPSCore's TYPED_TEST_CASE, which GoogleTest 1.10 deprecated.
    target_compile_options(${name} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
    alps_add_mpi_gtest(${name} COMPONENT utilities RANKS 2 OWN_MAIN)
    alps_add_mpi_gtest(${name} COMPONENT utilities RANKS 3 OWN_MAIN)
  endforeach()
endif()
