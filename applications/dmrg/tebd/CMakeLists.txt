#  Copyright Matthias Troyer 2009 - 2010.
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

find_package(HDF5)
if(LAPACK_FOUND AND ALPS_BUILD_FORTRAN)
  set_source_files_properties(
    tebd.f90
    PROPERTIES Fortran_PREPROCESS ON
  )
  set(SOURCES core/GlobalData.f90 core/LinearOps.f90 core/HamiOps.f90 
              core/StateOps.f90 core/LocalOps.f90 core/ObOps.f90
              core/PropOps.f90 core/Hdf5Interface.f90 tebd.f90)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(tebd ${SOURCES})
  target_link_libraries(tebd ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES} ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  if(HDF5_DEFINITIONS)
    foreach (_hdf5_def ${HDF5_DEFINITIONS})
      string(REPLACE "-D" "" _escaped_def ${_hdf5_def})
      list(APPEND ALPS_TEBD_DEFINITIONS ${_escaped_def})
    endforeach(_hdf5_def)
    set_target_properties(tebd PROPERTIES COMPILE_DEFINITIONS "${ALPS_TEBD_DEFINITIONS}")
  endif(HDF5_DEFINITIONS)
  if(ALPS_PYTHON_WHEEL)
    install(TARGETS tebd RUNTIME DESTINATION pyalps/bin COMPONENT applications)
    if(APPLE)
      set_target_properties(tebd PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
    else(APPLE)
      set_target_properties(tebd PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
    endif(APPLE)
  else()
    install(TARGETS tebd RUNTIME DESTINATION bin COMPONENT applications)
  endif()
else(LAPACK_FOUND AND ALPS_BUILD_FORTRAN)
  message(STATUS "tebd will not be built")
endif(LAPACK_FOUND AND ALPS_BUILD_FORTRAN)
