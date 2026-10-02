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
target_link_libraries(alps_hdf5 PUBLIC alps_headers ALPS::utilities
  HDF5::HDF5 Boost::filesystem Boost::thread Threads::Threads)
# FindHDF5 may expose provider targets or absolute paths through its aggregate.
get_target_property(hdf5_dependencies HDF5::HDF5 INTERFACE_LINK_LIBRARIES)
alps_dependency_rpath(alps_hdf5 ${hdf5_dependencies})
alps_copy_runtime_dependencies(alps_hdf5)
generate_export_header(alps_hdf5 EXPORT_MACRO_NAME ALPS_HDF5_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/hdf5_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_HDF5_STATIC_DEFINE)
endif()
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/hdf5_export.h")
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

install(TARGETS alps_hdf5 EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
