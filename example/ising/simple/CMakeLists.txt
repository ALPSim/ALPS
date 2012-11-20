#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_EXAMPLES)

    add_executable(example_ising_single single)
    add_dependencies(example_ising_single alps)
    target_link_libraries(example_ising_single alps)

endif (ALPS_BUILD_EXAMPLES)
