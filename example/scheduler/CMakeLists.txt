#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#include_directories(${PROJECT_BINARY_DIR}/src)
#include_directories(${PROJECT_SOURCE_DIR}/src)
#include_directories(${PROJECT_SOURCE_DIR}/tool)
#include_directories(${Boost_ROOT_DIR})

if (ALPS_BUILD_EXAMPLES)
  add_executable(ising1 main.C ising.C)
  add_executable(ising2 main2.C ising2.C)
  add_executable(ising3 main3.C ising.C)
  add_dependencies(ising1 alps)
  add_dependencies(ising2 alps)
  add_dependencies(ising3 alps)
  target_link_libraries(ising1 alps)
  target_link_libraries(ising2 alps)
  target_link_libraries(ising3 alps)
  #ADD_TEST(ising1 ising1 ./parm.in.xml)
  #ADD_TEST(ising2 ising2 ./parm2.in.xml)
  #ADD_TEST(ising3 ising3 ./parm3.in.xml)
endif (ALPS_BUILD_EXAMPLES)
