#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
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

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${Boost_ROOT_DIR})

IF(NOT ALPS_LLVM_WORKAROUND)
  FOREACH (name clone_info clone_phase clone_timer exmc_optimize exp_number footprint info_test integer_range linear_regression merge percentage temperature_scan time version wl_weight)
    add_executable(${name} ${name}.C)
    add_dependencies(${name} alps)
    target_link_libraries(${name} alps)
    add_alps_test(${name})
    set_property(TEST ${name} PROPERTY LABELS parapack)
  ENDFOREACH(name)

  IF(ALPS_HAVE_MPI)
    FOREACH (name collect_mpi comm_mpi filelock_mpi halt_mpi info_test_mpi process_mpi)
      add_executable(${name} ${name}.C)
      add_dependencies(${name} alps)
      target_link_libraries(${name} alps)
      # add_alps_test(${name})
      # set_property(TEST ${name} PROPERTY LABELS parapack)
    ENDFOREACH(name)
  ENDIF(ALPS_HAVE_MPI)
ENDIF(NOT ALPS_LLVM_WORKAROUND)
