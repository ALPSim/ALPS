#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
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

foreach(name os sizeof)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::osiris)
  alps_add_gtest(${name} COMPONENT osiris LABELS compatibility)
endforeach()

target_link_libraries(os PRIVATE ALPS::utilities)

add_executable(osiris_comm_contract comm_contract.cpp)
target_link_libraries(osiris_comm_contract PRIVATE ALPS::osiris)
alps_add_gtest(osiris_comm_contract COMPONENT osiris LABELS compatibility)

add_executable(osiris_dump_serialization dump_serialization.cpp)
target_link_libraries(osiris_dump_serialization PRIVATE ALPS::osiris)
# The Boost archive adapter in <alps/osiris/archivedump.h> is tested when
# Boost.Serialization is installed; ALPS itself does not need it.
find_package(Boost ${Boost_VERSION} EXACT CONFIG QUIET COMPONENTS serialization)
if(TARGET Boost::serialization)
  target_link_libraries(osiris_dump_serialization PRIVATE Boost::serialization)
  target_compile_definitions(osiris_dump_serialization PRIVATE ALPS_TEST_BOOST_ARCHIVES)
else()
  message(STATUS "Boost.Serialization not found: skipping the archivedump.h test")
endif()
alps_add_gtest(osiris_dump_serialization COMPONENT osiris LABELS compatibility)
