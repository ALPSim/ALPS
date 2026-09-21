cmake_minimum_required(VERSION 4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating single core ising program
add_executable(single ising.cpp single.cpp)
target_link_libraries(single PRIVATE ALPS::alps)

if(ALPS_MPI_ENABLED)

  # rule for generating mpi ising program
  add_executable(mpi ising.cpp mpi.cpp)
  target_link_libraries(mpi PRIVATE ALPS::alps)

  # rule for generating mpi pscan ising program
  add_executable(mpi_pscan ising.cpp mpi_pscan.cpp)
  target_link_libraries(mpi_pscan PRIVATE ALPS::alps)

endif(ALPS_MPI_ENABLED)
