#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

IF(ALPS_NGS_INCLUDE_NEW_ALEA_TEST_AND_EXAMPLE)
    FOREACH (name   example_accumulator 
                    example_accumulator_set
                    example_custom_accum 
                    example_histogram
                    example_vector_operators
                    )
        add_executable(${name} ${name}.cpp)
        target_link_libraries(${name} alps)
        add_alps_test(${name})
    ENDFOREACH(name)
ENDIF(ALPS_NGS_INCLUDE_NEW_ALEA_TEST_AND_EXAMPLE)
