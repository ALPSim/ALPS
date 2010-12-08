#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LPSolve_FOUND AND NOT MSVC)
  include_directories(${LPSolve_INCLUDE_DIR})
  add_executable(dirloop_sse_v1
                 SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp
		 SSE.Update.cpp SSE.cpp)
  target_link_libraries(dirloop_sse_v1 alps ${LPSolve_LIBRARIES})
  install(TARGETS dirloop_sse_v1 RUNTIME DESTINATION bin COMPONENT applications)
endif(LPSolve_FOUND AND NOT MSVC)
if (LPSolve_FOUND AND APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -idirafter /usr/include/malloc")
endif(LPSolve_FOUND AND APPLE)
add_executable(dirloop_sse_evaluate evaluate.C)
target_link_libraries(dirloop_sse_evaluate alps)
install(TARGETS dirloop_sse_evaluate RUNTIME DESTINATION bin COMPONENT applications)
