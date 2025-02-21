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


set(ALPS_INSPECT_SOURCES end_check.cpp license_check.cpp path_name_check.cpp 
          inspect.cpp tab_check.cpp)
set(INSPECT_SOURCES deprecated_macro_check.cpp link_check.cpp crlf_check.cpp  unnamed_namespace_check.cpp
          ascii_check.cpp copyright_check.cpp minmax_check.cpp apple_macro_check.cpp 
          assert_macro_check.cpp)
set(DIRECTORY "${Boost_ROOT_DIR}/tools/inspect")
include_directories(${DIRECTORY})
foreach(S ${INSPECT_SOURCES})
  if(EXISTS ${DIRECTORY}/${S})
    set(ALPS_INSPECT_SOURCES ${ALPS_INSPECT_SOURCES} ${DIRECTORY}/${S})
  endif(EXISTS ${DIRECTORY}/${S})
endforeach(S)

add_executable(alps_inspect ${ALPS_INSPECT_SOURCES})
target_link_libraries(alps_inspect boost)
install(TARGETS alps_inspect RUNTIME DESTINATION bin  COMPONENT tools)
