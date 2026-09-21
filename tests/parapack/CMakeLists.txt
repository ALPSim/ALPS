#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(
  name
  clone_info
  clone_phase
  clone_timer
  exmc_optimize
  exp_number
  footprint
  info_test
  integer_range
  linear_regression
  merge
  percentage
  temperature_scan
  time
  version
  wl_weight)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS parapack)
endforeach()

if(ALPS_ENABLE_MPI)
  foreach(name collect_mpi comm_mpi filelock_mpi halt_mpi info_test_mpi process_mpi)
    add_executable(${name} ${name}.C)

    target_link_libraries(${name} PRIVATE ALPS::alps)
    # alps_add_test(${name})
    # set_property(TEST ${name} PROPERTY LABELS parapack)
  endforeach()
endif()
