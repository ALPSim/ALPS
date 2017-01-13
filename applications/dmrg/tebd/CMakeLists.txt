#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.f90rg/LICENSE_1_0.txt)

find_package(HDF5)
if(LAPACK_FOUND AND ALPS_BUILD_FORTRAN AND HDF5_STATIC_LIBRARIES)
  set(SOURCES core/GlobalData.f90 core/LinearOps.f90 core/HamiOps.f90 
              core/StateOps.f90 core/LocalOps.f90 core/ObOps.f90
              core/PropOps.f90 core/Hdf5Interface.f90 tebd.f90)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(tebd ${SOURCES})
  target_link_libraries(tebd ${HDF5_STATIC_LIBRARIES} ${HDF5_LIBRARIES} ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  if(HDF5_DEFINITIONS)
    foreach (_hdf5_def ${HDF5_DEFINITIONS})
      string(REPLACE "-D" "" _escaped_def ${_hdf5_def})
      list(APPEND ALPS_TEBD_DEFINITIONS ${_escaped_def})
    endforeach(_hdf5_def)
    set_target_properties(tebd PROPERTIES COMPILE_DEFINITIONS "${ALPS_TEBD_DEFINITIONS}")
  endif(HDF5_DEFINITIONS)
  install(TARGETS tebd RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND AND ALPS_BUILD_FORTRAN)
  message(STATUS "tebd will not be built")
endif(LAPACK_FOUND AND ALPS_BUILD_FORTRAN AND HDF5_STATIC_LIBRARIES)
