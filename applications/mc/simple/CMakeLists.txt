#  Copyright Synge Todo 2010 - 2015.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_executable(simplemc main.C evaluator.C ising.C xy.C heisenberg.C)
target_link_libraries(simplemc alps)
install(TARGETS simplemc RUNTIME DESTINATION bin COMPONENT applications)

enable_testing()
add_alps_test(simplemc_ising simplemc ising ising)
add_alps_test(simplemc_xy simplemc xy xy)
add_alps_test(simplemc_heisenberg simplemc heisenberg heisenberg)
