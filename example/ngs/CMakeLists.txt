#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)

    foreach(program params)
        add_executable(ngs_${program} ${program}.cpp)
        target_link_libraries(ngs_${program} alps)
    endforeach(program)

endif (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)
