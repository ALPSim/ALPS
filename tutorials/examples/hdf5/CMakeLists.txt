#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

#
# examples using the hdf5 library
#

foreach(program enum_as_class enum_vectorizable pair_int_vectorizable)
  add_executable(${program} ${program}.cpp)

  target_link_libraries(${program} PRIVATE ALPS::alps)
  alps_add_test(${program})
endforeach()
