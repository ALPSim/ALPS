#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
#   SPDX-License-Identifier: MIT

foreach(name ngs_hdf5 ngs_hash)
  add_executable(${name} ${name}.cpp)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
endforeach()

#TODO: does not compile in clang
add_subdirectory(params)
add_subdirectory(scheduler)
