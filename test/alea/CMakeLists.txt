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

foreach(
  name
  complexobservable
  detailedbinning
  histogram
  histogram2
  signed
  simpleobseval
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
