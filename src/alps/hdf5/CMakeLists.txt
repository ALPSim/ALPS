# Module-owned implementation and public headers; linked into the common ALPS runtime.
target_sources(alps PRIVATE
  src/archive.cpp
  src/signal.cpp)
target_sources(alps_headers INTERFACE
  FILE_SET hdf5_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/hdf5/archive.hpp
    include/alps/hdf5/array.hpp
    include/alps/hdf5/complex.hpp
    include/alps/hdf5/errors.hpp
    include/alps/hdf5/map.hpp
    include/alps/hdf5/matrix.hpp
    include/alps/hdf5/multi_array.hpp
    include/alps/hdf5/numeric_vector.hpp
    include/alps/hdf5/pair.hpp
    include/alps/hdf5/pointer.hpp
    include/alps/hdf5/shared_array.hpp
    include/alps/hdf5/stdarray.hpp
    include/alps/hdf5/tuple.hpp
    include/alps/hdf5/ublas/matrix.hpp
    include/alps/hdf5/ublas/vector.hpp
    include/alps/hdf5/valarray.hpp
    include/alps/hdf5/vector.hpp
    include/alps/hdf5.hpp
    include/alps/ngs/signal.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
