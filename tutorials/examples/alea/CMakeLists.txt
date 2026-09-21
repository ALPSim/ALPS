#  Copyright Matthias Troyer and Maximilian Poprawe 2011 - 2012.
#   SPDX-License-Identifier: MIT

configure_file(testfile.h5 testfile.h5 COPYONLY)

foreach(example mean variance autocorrelation error running_mean)
  add_executable(example_${example} example_${example}.cpp)

  target_link_libraries(example_${example} PRIVATE ALPS::alps)
endforeach()
