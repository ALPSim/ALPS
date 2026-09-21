#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(dirloop_sse main.cc lp_sse.cpp)
target_link_libraries(dirloop_sse PRIVATE ALPS::alps)
alps_install_applications(dirloop_sse)
