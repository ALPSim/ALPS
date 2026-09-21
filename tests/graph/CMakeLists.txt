#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

if(ALPS_BUILD_EXTENSIVE_TESTS)

  foreach(
    name
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
    colored_lattice_constant_test2)
    add_executable(${name} ${name}.cpp)

    target_link_libraries(${name} PRIVATE ALPS::alps)
    alps_add_test(${name})
    set_property(TEST ${name} PROPERTY LABELS "graph;extensive")
  endforeach()
endif()
