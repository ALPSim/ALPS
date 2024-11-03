#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
if(ALPS_LLVM_WORKAROUND)
  message(STATUS "dmrg will not be built since llvm workaround is enabled")
else(ALPS_LLVM_WORKAROUND)
  add_definitions(${LAPACK_DEFINITIONS} -DBIND_FORTRAN_LOWERCASE_UNDERSCORE)
  include_directories(${PROJECT_SOURCE_DIR}/applications/dmrg/dmrg)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(dmrg dmrg.C factory.C)
  target_link_libraries(dmrg alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
if(ALPS_PYTHON_WHEEL)
  install(TARGETS dmrg RUNTIME DESTINATION pyalps/bin COMPONENT applications)
set_target_properties(dmrg  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
else()
  install(TARGETS dmrg RUNTIME DESTINATION bin COMPONENT applications)
endif()
endif(ALPS_LLVM_WORKAROUND)
else(LAPACK_FOUND)
  message(STATUS "dmrg will not be built since lapack library is not found")
endif(LAPACK_FOUND)
