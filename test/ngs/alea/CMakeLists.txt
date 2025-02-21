#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

IF(ALPS_NGS_USE_NEW_ALEA)
    FOREACH (name ngs_alea_next ngs_alea_compare)
        add_executable(${name} ${name}.cpp)
        add_dependencies(${name} alps)
        target_link_libraries(${name} alps)
        add_alps_test(${name})
    ENDFOREACH(name)
ENDIF(ALPS_NGS_USE_NEW_ALEA)

IF(FALSE AND ALPS_NGS_USE_NEW_ALEA)
    FOREACH (name   ngs_alea_count_test_compile
                    ngs_alea_count_test_runtime
                    ngs_alea_ctor_test_compile
                    ngs_alea_ctor_test_runtime
                    ngs_alea_error_test_compile
                    ngs_alea_error_test_runtime
                    ngs_alea_fix_size_test_compile
                    ngs_alea_fix_size_test_runtime
                    ngs_alea_max_num_test_compile
                    ngs_alea_max_num_test_runtime
                    ngs_alea_log_test_compile
                    ngs_alea_log_test_runtime
                    ngs_alea_mean_test_runtime
                    ngs_alea_mean_test_compile
                    ngs_alea_stream_test_compile
                    ngs_alea_stream_test_runtime
                    ngs_alea_value_type_test
                    ngs_alea_weight_type_test
                    ngs_alea_wrapper_test_compile
                    ngs_alea_wrapper_test_runtime
    )
        add_executable(${name} ${name}.cpp)
        add_dependencies(${name} alps)
        target_link_libraries(${name} alps)
        add_alps_test(${name})
    ENDFOREACH(name)
ENDIF(FALSE AND ALPS_NGS_USE_NEW_ALEA)
