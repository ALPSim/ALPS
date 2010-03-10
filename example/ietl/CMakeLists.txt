#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# examples of using the IETL
# 

if (ALPS_BUILD_EXAMPLES)
  set (EXAMPLES arnoldi1)
     
  foreach(program ${EXAMPLES})
    add_executable(${program} ${program}.C)
    target_link_libraries(${program} alps)
  endforeach(program)
endif (ALPS_BUILD_EXAMPLES)

