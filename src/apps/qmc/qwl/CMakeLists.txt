#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(qwl qwl.C)
add_executable(qwl_evaluate qwl_evaluate.C)

target_link_libraries(qwl PRIVATE ALPS::alps)
target_link_libraries(qwl_evaluate PRIVATE ALPS::alps)
alps_install_applications(qwl qwl_evaluate)
