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
include_directories(${Boost_ROOT_DIR})

IF(NOT ALPS_LLVM_WORKAROUND)
  FOREACH (name complexobservable detailedbinning histogram histogram2 signed simpleobseval testobservableset observableset_hdf5 observableset_xml vectorobseval mcdata mcdata2 mcanalyze mergeobs_plus mergeobs_plus2 mergeobs_times mergeobs_times2 mergeobs_divide mergeobs_log)
    add_executable(${name} ${name}.C)
    add_dependencies(${name} alps)
    target_link_libraries(${name} alps)
    add_alps_test(${name})
  ENDFOREACH(name)
ENDIF(NOT ALPS_LLVM_WORKAROUND)

add_executable(dumpbench dumpbench.C)
target_link_libraries(dumpbench alps ${HDF5_LIBRARIES})
if(HDF5_DEFINITIONS)
   foreach (_hdf5_def ${HDF5_DEFINITIONS})
     string(REPLACE "-D" "" _escaped_def ${_hdf5_def})
     list(APPEND ALPS_TEBD_DEFINITIONS ${_escaped_def})
   endforeach(_hdf5_def)
   set_target_properties(dumpbench PROPERTIES COMPILE_DEFINITIONS "${ALPS_TEBD_DEFINITIONS}")
endif(HDF5_DEFINITIONS)