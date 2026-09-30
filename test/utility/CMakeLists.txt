#  Copyright Synge Todo 2012.
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

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

FOREACH(name vmusage bitops)
  add_executable(${name} ${name}.cpp)
  add_dependencies(${name} alps)
  target_link_libraries(${name} alps)
  add_alps_test(${name})
  set_property(TEST ${name} PROPERTY LABELS utility)
ENDFOREACH(name)

add_executable(dmrg_temporary_files dmrg_temporary_files.cpp)
target_include_directories(dmrg_temporary_files PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(dmrg_temporary_files alps)
add_test(NAME dmrg_temporary_files COMMAND dmrg_temporary_files)
set_tests_properties(dmrg_temporary_files PROPERTIES LABELS "utility;dmrg")

if(TARGET dmrg AND TARGET parameter2xml)
  add_test(NAME dmrg_scratch_cleanup
    COMMAND ${CMAKE_COMMAND}
      -Ddmrg=$<TARGET_FILE:dmrg>
      -Dparameter2xml=$<TARGET_FILE:parameter2xml>
      -Dsource_dir=${PROJECT_SOURCE_DIR}
      -Dxml_dir=${PROJECT_BINARY_DIR}/lib/xml
      -Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/dmrg-scratch
      -P ${CMAKE_CURRENT_SOURCE_DIR}/dmrg_scratch_cleanup.cmake)
  if(UNIX)
    add_test(NAME dmrg_scratch_cleanup_long_paths
      COMMAND ${CMAKE_COMMAND}
        -Ddmrg=$<TARGET_FILE:dmrg>
        -Dparameter2xml=$<TARGET_FILE:parameter2xml>
        -Dsource_dir=${PROJECT_SOURCE_DIR}
        -Dxml_dir=${PROJECT_BINARY_DIR}/lib/xml
        -Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/dmrg-long-scratch
        -Dlong_relative_paths=ON
        -P ${CMAKE_CURRENT_SOURCE_DIR}/dmrg_scratch_cleanup.cmake)
    set_tests_properties(dmrg_scratch_cleanup_long_paths PROPERTIES LABELS "dmrg" TIMEOUT 120)
  endif()
  set_tests_properties(dmrg_scratch_cleanup PROPERTIES LABELS "dmrg" TIMEOUT 120)
endif()
