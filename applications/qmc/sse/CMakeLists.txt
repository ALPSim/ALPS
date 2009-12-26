if(LPSolve_FOUND)
  include_directories(${LPSolve_INCLUDE_DIR})
  add_executable(dirloop_sse
                 SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp
		 SSE.Update.cpp SSE.cpp)
  target_link_libraries(dirloop_sse alps boost ${LPSolve_LIBRARIES})
  install(TARGETS dirloop_sse RUNTIME DESTINATION bin)
else(LPSolve_FOUND)
  message(STATUS "sse will not be built since lp_sovle library is not found")
endif(LPSolve_FOUND)
