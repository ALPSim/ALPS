#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

# Looper headers use their package-relative include paths.
add_executable(loop loop.C loop_custom.C loop_model.C path_integral.C sse.C)
target_include_directories(loop PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(loop PRIVATE ALPS::alps)
alps_install_applications(loop)
