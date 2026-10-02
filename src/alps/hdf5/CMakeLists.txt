# Archive state and its signal cleanup handler share one runtime owner.
add_library(alps_hdf5
  src/archive.cpp
  src/signal.cpp)
add_library(ALPS::hdf5 ALIAS alps_hdf5)
set_target_properties(alps_hdf5 PROPERTIES
  EXPORT_NAME hdf5
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_hdf5 PUBLIC ALPS::configuration ALPS::containers ALPS::utilities
  HDF5::HDF5 Boost::filesystem Boost::thread Threads::Threads)
# FindHDF5 may expose provider targets or absolute paths through its aggregate.
get_target_property(hdf5_dependencies HDF5::HDF5 INTERFACE_LINK_LIBRARIES)
alps_dependency_rpath(alps_hdf5 ${hdf5_dependencies})
alps_copy_runtime_dependencies(alps_hdf5)
generate_export_header(alps_hdf5 EXPORT_MACRO_NAME ALPS_HDF5_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/hdf5_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_hdf5 PUBLIC ALPS_HDF5_STATIC_DEFINE)
endif()
target_sources(alps_hdf5 PUBLIC
  FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include"
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/hdf5_export.h")
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

alps_register_module(NAME hdf5 HEADER_SET hdf5_headers
  HEADER_TARGETS alps_hdf5
  TARGETS alps alps_hdf5
  DEPENDS configuration containers utilities
  GENERATED_HEADERS alps/hdf5_export.h
)
