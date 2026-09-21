#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(dmrg dmrg.C factory.C)
target_include_directories(dmrg PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(dmrg PRIVATE ALPS::alps)
alps_install_applications(dmrg)
