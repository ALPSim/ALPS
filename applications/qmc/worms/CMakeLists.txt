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
 
add_library(worm_impl STATIC WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C)
target_link_libraries(worm_impl alps)

add_executable(worm main.C)
add_executable(worm_evaluate evaluate.C)

target_link_libraries(worm worm_impl)
target_link_libraries(worm_evaluate worm_impl)
if(ALPS_PYTHON_WHEEL)
  if(APPLE)
    set_target_properties(worm PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
    set_target_properties(worm_evaluate PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
  else(APPLE)
    set_target_properties(worm PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
    set_target_properties(worm_evaluate PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
  endif(APPLE)
  install(TARGETS worm RUNTIME DESTINATION pyalps/bin COMPONENT applications)
  install(TARGETS worm_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
  install(TARGETS worm RUNTIME DESTINATION bin COMPONENT applications)
  install(TARGETS worm_evaluate RUNTIME DESTINATION bin COMPONENT applications)
endif()