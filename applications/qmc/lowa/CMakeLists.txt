#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(lowa main.cpp lowa.cpp)
  target_link_libraries(lowa alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  install(TARGETS lowa RUNTIME DESTINATION bin COMPONENT applications)  
  if (NOT WIN32)

    add_executable(extract_timeseries_hdf5_measurements extract_timeseries_hdf5_measurements.cpp)
    add_executable(merge_broken_timeseries_hdf5_measurements merge_broken_timeseries_hdf5_measurements.cpp)  
    add_executable(from_3d_to_2d_hdf5_measurements from_3d_to_2d_hdf5_measurements.cpp)
    add_executable(from_2d_to_1d_cylindrically_binned_hdf5_measurements from_2d_to_1d_cylindrically_binned_hdf5_measurements.cpp)
    add_executable(evaluate_observable_observable_xi_coupling_timeseries_2d_hdf5_measurements evaluate_observable_observable_xi_coupling_timeseries_2d_hdf5_measurements.cpp)

    target_link_libraries(extract_timeseries_hdf5_measurements alps)
    target_link_libraries(merge_broken_timeseries_hdf5_measurements alps)
    target_link_libraries(from_3d_to_2d_hdf5_measurements alps)
    target_link_libraries(from_2d_to_1d_cylindrically_binned_hdf5_measurements alps)
    target_link_libraries(evaluate_observable_observable_xi_coupling_timeseries_2d_hdf5_measurements alps)

    install(TARGETS extract_timeseries_hdf5_measurements  merge_broken_timeseries_hdf5_measurements from_3d_to_2d_hdf5_measurements from_2d_to_1d_cylindrically_binned_hdf5_measurements evaluate_observable_observable_xi_coupling_timeseries_2d_hdf5_measurements RUNTIME DESTINATION bin COMPONENT applications)  
  endif(NOT WIN32)
else(LAPACK_FOUND)
  message(STATUS "lowa will not be built since lapack library has not been found")
endif(LAPACK_FOUND)
