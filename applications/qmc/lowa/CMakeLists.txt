#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(lowa main.cpp lowa.cpp)
  add_executable(extract_timeseries_hdf5_measurements extract_timeseries_hdf5_measurements.cpp)
  add_executable(merge_broken_timeseries_hdf5_measurements merge_broken_timeseries_hdf5_measurements.cpp)  
  target_link_libraries(lowa alps boost  ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  target_link_libraries(extract_timeseries_hdf5_measurements alps boost)
  target_link_libraries(merge_broken_timeseries_hdf5_measurements alps boost)
  install(TARGETS lowa RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS extract_timeseries_hdf5_measurements  RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS merge_broken_timeseries_hdf5_measurements RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND)
  message(STATUS "lowa will not be built since lapack library has not been found")
endif(LAPACK_FOUND)
