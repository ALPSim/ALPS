#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


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
