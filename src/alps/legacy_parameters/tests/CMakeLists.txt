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

foreach(name expression expression2 flatten parameter parameterlist parameters)
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_gtest(${name} COMPONENT legacy_parameters
    PROPERTIES ENVIRONMENT "DIR=/home/alps")
endforeach()

# XML serialization spelling and ordering remain explicit format contracts.
foreach(name parameterlist_xml parameters_xml)
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_gtest(${name} COMPONENT legacy_parameters LABELS compatibility format)
endforeach()

add_executable(parameters_hdf5 parameters_hdf5.C)
target_link_libraries(parameters_hdf5 PRIVATE ALPS::alps)
alps_add_gtest(parameters_hdf5 COMPONENT legacy_parameters LABELS integration hdf5)

if(ALPS_ENABLE_MPI)
  add_executable(parameters_mpi parameters_mpi.C)
  target_link_libraries(parameters_mpi PRIVATE ALPS::alps)
  alps_add_mpi_gtest(parameters_mpi COMPONENT legacy_parameters RANKS 2)
  alps_add_mpi_gtest(parameters_mpi COMPONENT legacy_parameters RANKS 3)
endif()
