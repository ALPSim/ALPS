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
else(LAPACK_FOUND)
  message(STATUS "lowa will not be built since lapack library has not been found")
endif(LAPACK_FOUND)
