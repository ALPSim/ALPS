#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_library(worm_impl STATIC WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C)
target_link_libraries(worm_impl PUBLIC ALPS::alps)

add_executable(worm main.C)
add_executable(worm_evaluate evaluate.C)

target_link_libraries(worm PRIVATE worm_impl)
target_link_libraries(worm_evaluate PRIVATE worm_impl)
alps_install_applications(worm worm_evaluate)
