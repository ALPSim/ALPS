#  Copyright Matthias Troyer and Maximilian Poprawe 2011 - 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_SOURCE_DIR}/tool/alea)
#include_directories(${PROJECT_BINARY_DIR}/src)
#include_directories(${PROJECT_SOURCE_DIR}/src)
#include_directories(${Boost_ROOT_DIR})


  FOREACH(program mean variance)
    add_executable(${program} ${program}.cpp)
    add_dependencies(${program} alps)
    target_link_libraries(${program} alps)
    #install(TARGETS ${program} RUNTIME DESTINATION bin COMPONENT tools) 
  ENDFOREACH(program)


