#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# examples using the hdf5 library
# 

if (ALPS_BUILD_EXAMPLES)
    foreach(program enum_as_class enum_vectorizable pair_int_vectorizable)
      add_executable(${program} ${program}.cpp)
      add_dependencies(${program} alps)
      target_link_libraries(${program} alps)
      add_alps_test(${program})
    endforeach(program)
endif (ALPS_BUILD_EXAMPLES)

