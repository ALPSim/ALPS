#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_executable(checksign checksign.C)
target_link_libraries(checksign alps)
install(TARGETS checksign RUNTIME DESTINATION bin COMPONENT applications)
