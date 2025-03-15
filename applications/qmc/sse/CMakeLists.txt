#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

if(LPSolve_FOUND AND NOT MSVC)
  include_directories(${LPSolve_INCLUDE_DIR})
  add_executable(dirloop_sse_v1
                 SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp
		 SSE.Update.cpp SSE.cpp)
  target_link_libraries(dirloop_sse_v1 alps ${LPSolve_LIBRARIES})
  if(ALPS_PYTHON_WHEEL)
  set_target_properties(dirloop_sse_v1  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
  install(TARGETS dirloop_sse_v1 RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  else()
  install(TARGETS dirloop_sse_v1 RUNTIME DESTINATION bin COMPONENT applications)
  endif()
endif(LPSolve_FOUND AND NOT MSVC)
if (LPSolve_FOUND AND APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -idirafter /usr/include/malloc")
endif(LPSolve_FOUND AND APPLE)
add_executable(dirloop_sse_evaluate evaluate.C)
target_link_libraries(dirloop_sse_evaluate alps)
if(ALPS_PYTHON_WHEEL)
set_target_properties(dirloop_sse_evaluate  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
install(TARGETS dirloop_sse_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS dirloop_sse_evaluate RUNTIME DESTINATION bin COMPONENT applications)
endif()