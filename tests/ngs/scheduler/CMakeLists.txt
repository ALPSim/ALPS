#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
#   SPDX-License-Identifier: MIT

if(ALPS_ENABLE_MPI)

  add_executable(scheduler_sum_mpi sum_mpi.cpp)

  target_link_libraries(scheduler_sum_mpi PRIVATE ALPS::alps)

endif(ALPS_ENABLE_MPI)

add_executable(scheduler_sum_single sum_single.cpp)

target_link_libraries(scheduler_sum_single PRIVATE ALPS::alps)
