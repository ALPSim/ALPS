#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(
  name
  complexobservable
  detailedbinning
  histogram
  histogram2
  signed
  simpleobseval
  simpleobs_initialization
  testobservableset
  observableset_hdf5
  observableset_xml
  vectorobseval
  mcdata
  mcdata2
  mcdata_transform_variance
  mcanalyze
  mergeobs_plus
  mergeobs_plus2
  mergeobs_times
  mergeobs_times2
  mergeobs_divide
  mergeobs_log)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
endforeach()

add_executable(dumpbench dumpbench.C)
target_link_libraries(dumpbench PRIVATE ALPS::alps)
