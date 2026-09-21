#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(checksign checksign.C)
target_link_libraries(checksign PRIVATE ALPS::alps)
alps_install_applications(checksign)
