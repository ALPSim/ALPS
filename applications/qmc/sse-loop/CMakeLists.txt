#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

 add_executable(sse-loop loop.C sse.C sse.h vertex.C vertex.h)
 target_link_libraries(sse-loop alps boost)
