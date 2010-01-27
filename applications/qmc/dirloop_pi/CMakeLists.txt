#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_executable(dirloop_pi mail.cpp lowa.cpp)
target_link_libraries(dirloop_pi alps boost)
install(TARGETS dirloop_pi RUNTIME DESTINATION bin COMPONENT applications)
