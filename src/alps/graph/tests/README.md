# Graph test migration inventory

All 17 formerly registered extensive graph executables now contain GoogleTest
assertions and register individual cases with CTest (`graph.<target>.*`). Their
input graphs, enumerated integer references, random seed 23, 50 random graphs,
and 100 isomorphic relabelings per random graph are preserved.

| Former test | Current assertion contract |
| --- | --- |
| `canonical_label_test` | Three cases check each historical canonical label encoding and the expected equality/inequality relations. |
| `canonical_label_with_color_symmetries_test` | Seven cases check label encodings, symmetry relations, and explicit color maps. |
| `canonical_label_random_graphs_test` | Every relabeling preserves its canonical label, with graph/relabeling indices in failure diagnostics. |
| `orbit_test` | Five cases compare each orbit's actual vertex membership with the historical reference partition. |
| `is_embeddable_with_color_symmetries_test` | Four cases assert every previously printed embeddability boolean, including pinned vertices. |
| `utils_test` | All 720 color mappings exist, are unique, and respect their color partitions. |
| `iso_simple` | Existing hand-constructed examples check canonical-order permutations, known isomorphisms and non-isomorphisms, and nonempty subgraph generation. This was previously print-only smoke coverage. |
| `lattice_constant_square_test`, `lattice_constant_tri_test` | Every independent integer lattice-constant reference is asserted. |
| `colored_lattice_constant_test`, `colored_lattice_constant_test2` | Every colored-graph lattice-constant reference is asserted. |
| `lattice_constant_matrix` | Every output vertex is compared with the preexisting reference calculation. |
| `subgraph_generator_test`, `subgraph_generator_test_colored_edges`, `subgraph_generator_test_colored_edges_with_sym` | Integer counts at every edge order and their totals come from the existing fixtures. |
| `subgraph_generator_test_colored_edges2` | Generated graphs have allowed edge counts and the result is nonempty; this previously had no output reference. |
| `subgraph_generator_test_colored_edges_with_sym2` | Canonical graph sets agree with and without the color-symmetry optimization for both lattice constructions. |

The old `.output` files for these tests were retired after their expectations
were transferred into assertions. Canonical label encodings remain exact string
comparisons because they encode discrete graph data; orbit membership, mappings,
booleans, lattice constants, and graph counts are compared as structured values.
`graph_assertions.hpp` contains only thin assertion helpers, not a test runner.

`embedding_test.cpp` and its historical fixture remain unregistered: they refer
to the removed `alps/graph/embedding.hpp` and `embedding_iterator` API. There is
no runtime coverage to migrate without selecting a replacement public counting
API. Current embeddability and lattice-constant tests are active; claiming the
obsolete finite-graph embedding counts are covered would be inaccurate.
