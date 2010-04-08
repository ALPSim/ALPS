#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LPSolve_FOUND)
  include_directories(${LPSolve_INCLUDE_DIR})
  add_executable(dirloop_sse
                 SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp
		 SSE.Update.cpp SSE.cpp)
  target_link_libraries(dirloop_sse alps ${LPSolve_LIBRARIES})
  install(TARGETS dirloop_sse RUNTIME DESTINATION bin COMPONENT applications)
else(LPSolve_FOUND)
  message(STATUS "dirloop_sse will not be built since the lp_solve library has not found")
endif(LPSolve_FOUND)
