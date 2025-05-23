#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
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

option(ALPS_BUILD_GRAPH_TESTS "Build the graph tests" OFF)
mark_as_advanced(ALPS_BUILD_GRAPH_TESTS)

IF (ALPS_BUILD_GRAPH_TESTS)

    FOREACH (name
            utils_test
            iso_simple
            canonical_label_test
            canonical_label_random_graphs_test
            canonical_label_with_color_symmetries_test
            orbit_test
            is_embeddable_with_color_symmetries_test
            subgraph_generator_test
            subgraph_generator_test_colored_edges
            subgraph_generator_test_colored_edges2
            subgraph_generator_test_colored_edges_with_sym
            subgraph_generator_test_colored_edges_with_sym2
            lattice_constant_square_test
            lattice_constant_tri_test
            lattice_constant_matrix
            colored_lattice_constant_test
            colored_lattice_constant_test2
    )
      add_executable(${name} ${name}.cpp)
      add_dependencies(${name} alps)
      target_link_libraries(${name} alps)
      add_alps_test(${name})
      set_property(TEST ${name} PROPERTY LABELS graph)
    ENDFOREACH(name)
ENDIF(ALPS_BUILD_GRAPH_TESTS)
