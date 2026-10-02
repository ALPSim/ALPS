# Lattice source ownership.
target_sources(alps PRIVATE
  src/lattice/disorder.C
  src/lattice/latticedescriptor.C
  src/lattice/latticegraphdescriptor.C
  src/lattice/latticelibrary.C
  src/lattice/unitcell.C)

target_sources(alps_headers INTERFACE
  FILE_SET lattice_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/lattice.h
    include/alps/lattice/bond_compare.h
    include/alps/lattice/boundary.h
    include/alps/lattice/cell_traits.h
    include/alps/lattice/coordinate_traits.h
    include/alps/lattice/coordinategraph.h
    include/alps/lattice/coordinatelattice.h
    include/alps/lattice/dimensional_traits.h
    include/alps/lattice/disorder.h
    include/alps/lattice/graph.h
    include/alps/lattice/graph_helper.h
    include/alps/lattice/graph_traits.h
    include/alps/lattice/graphproperties.h
    include/alps/lattice/hypercubic.h
    include/alps/lattice/lattice.h
    include/alps/lattice/latticedescriptor.h
    include/alps/lattice/latticegraph.h
    include/alps/lattice/latticegraphdescriptor.h
    include/alps/lattice/latticelibrary.h
    include/alps/lattice/parity.h
    include/alps/lattice/point_traits.h
    include/alps/lattice/propertymap.h
    include/alps/lattice/simplecell.h
    include/alps/lattice/simplelattice.h
    include/alps/lattice/unitcell.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME lattice HEADER_SET lattice_headers
  TARGETS alps
  DEPENDS configuration containers expression legacy_parameters numerics osiris random utilities xml
)
