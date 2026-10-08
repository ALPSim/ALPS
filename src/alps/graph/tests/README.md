# Graph tests

Extensive tests compare canonical labels, orbit membership, color mappings,
embeddability, lattice constants and subgraph counts with exact integer/string
references. Random relabeling uses seed 23, 50 graphs and 100 relabelings per graph.
Color-symmetry optimizations must preserve the generated canonical graph sets.

`embedding_test.cpp` remains unregistered because it refers to the removed
`embedding.hpp` API; it provides no coverage of current graph interfaces.
