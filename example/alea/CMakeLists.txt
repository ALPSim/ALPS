#  Copyright Matthias Troyer and Maximilian Poprawe 2011 - 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#include_directories(${PROJECT_SOURCE_DIR}/tool/alea)
#include_directories(${PROJECT_BINARY_DIR}/src)
#include_directories(${PROJECT_SOURCE_DIR}/src)
#include_directories(${Boost_ROOT_DIR})


  FILE ( COPY testfile.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


  FOREACH(example mean variance autocorrelation error running_mean)
    add_executable(example_${example} example_${example}.cpp)
    add_dependencies(example_${example} alps)
    target_link_libraries(example_${example} alps)
    #install(TARGETS example_${example} RUNTIME DESTINATION bin COMPONENT example) 
  ENDFOREACH(example)

  
