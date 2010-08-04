#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LPSolve_FOUND)
  include_directories(${LPSolve_INCLUDE_DIR})
  add_executable(sse4 main.cc)
  target_link_libraries(sse4 alps ${LPSolve_LIBRARIES})
  install(TARGETS sse4 RUNTIME DESTINATION bin COMPONENT applications)
else(LPSolve_FOUND)
  message(STATUS "sse4 will not be built since the lp_solve library has not found")
endif(LPSolve_FOUND)
