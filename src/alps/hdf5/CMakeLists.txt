# Archive state and its signal cleanup handler share one runtime owner.
add_library(alps_hdf5
  src/archive.cpp
  src/signal.cpp)
target_link_libraries(alps_hdf5 PUBLIC ALPS::configuration ALPS::containers ALPS::utilities
  HDF5::HDF5 Boost::filesystem Boost::thread Threads::Threads)
# FindHDF5 may expose provider targets or absolute paths through its aggregate.
get_target_property(hdf5_dependencies HDF5::HDF5 INTERFACE_LINK_LIBRARIES)
alps_configure_component(alps_hdf5 hdf5 ${hdf5_dependencies})
target_sources(alps_hdf5 PUBLIC
  FILE_SET hdf5_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/hdf5.hpp
    include/alps/hdf5/archive.hpp
    include/alps/hdf5/array.hpp
    include/alps/hdf5/complex.hpp
    include/alps/hdf5/errors.hpp
    include/alps/hdf5/map.hpp
    include/alps/hdf5/multi_array.hpp
    include/alps/hdf5/pair.hpp
    include/alps/hdf5/pointer.hpp
    include/alps/hdf5/shared_array.hpp
    include/alps/hdf5/stdarray.hpp
    include/alps/hdf5/tuple.hpp
    include/alps/hdf5/ublas/matrix.hpp
    include/alps/hdf5/ublas/vector.hpp
    include/alps/hdf5/valarray.hpp
    include/alps/hdf5/vector.hpp
    include/alps/ngs/signal.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
