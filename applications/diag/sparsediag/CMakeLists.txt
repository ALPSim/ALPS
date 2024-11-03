#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(sparsediag sparsediag.C factory.C)
  target_link_libraries(sparsediag alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  if(ALPS_PYTHON_WHEEL)
  set_target_properties(sparsediag  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
  install(TARGETS sparsediag RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  else()
  install(TARGETS sparsediag RUNTIME DESTINATION bin COMPONENT applications)
  endif()
else(LAPACK_FOUND)
  message(STATUS "sparsediag will not be built since lapack library is not found")
endif(LAPACK_FOUND)
