#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_library(spinmc_impl STATIC abstract_fitter.C base_incr_fitter.C dummy_fitter.C
                               est_grad_fitter.C fitting_scheduler.C spinmc_factory.C)
target_link_libraries(spinmc_impl PUBLIC ALPS::alps)
add_executable(spinmc spinmc.C)
add_executable(spinmc_evaluate spinmc_evaluate.C)
target_link_libraries(spinmc PRIVATE spinmc_impl)
target_link_libraries(spinmc_evaluate PRIVATE spinmc_impl)
alps_install_applications(spinmc spinmc_evaluate)
