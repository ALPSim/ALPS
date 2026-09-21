#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_library(fulldiag_impl STATIC factory.C)
target_link_libraries(fulldiag_impl PUBLIC ALPS::alps)

add_executable(fulldiag fulldiag.C)
add_executable(fulldiag_evaluate fulldiag_evaluate.C)
target_link_libraries(fulldiag PRIVATE fulldiag_impl)
target_link_libraries(fulldiag_evaluate PRIVATE fulldiag_impl)

if(BUILD_TESTING)
  add_test(NAME fulldiag_evaluate_help COMMAND fulldiag_evaluate --help)
  set_tests_properties(fulldiag_evaluate_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")
  add_test(NAME fulldiag_evaluate_license COMMAND fulldiag_evaluate --license)
  set_tests_properties(fulldiag_evaluate_license PROPERTIES PASS_REGULAR_EXPRESSION
                                                            "Licensed under the MIT License")
  add_test(NAME fulldiag_evaluate_no_arguments COMMAND fulldiag_evaluate)
  set_tests_properties(fulldiag_evaluate_no_arguments PROPERTIES WILL_FAIL TRUE)
endif()

alps_install_applications(fulldiag fulldiag_evaluate)
