# Checkpoint serialization and communication state have one runtime owner.
add_library(alps_osiris
  src/osiris/buffer.C
  src/osiris/comm.C
  src/osiris/dump.C
  src/osiris/dumparchive.C
  src/osiris/mpdump.C
  src/osiris/process.C
  src/osiris/xdrdump.C
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_array.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_float.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_stdio.c")
target_compile_features(alps_osiris PRIVATE c_std_11)
target_link_libraries(alps_osiris PUBLIC ALPS::configuration Boost::serialization Boost::filesystem)
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_osiris PUBLIC MPI::MPI_CXX)
endif()
alps_configure_component(alps_osiris osiris)

target_sources(alps_osiris PUBLIC
  FILE_SET osiris_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/third_party/xdr/include"
  FILES
    include/alps/osiris.h
    include/alps/osiris/archivedump.h
    include/alps/osiris/boost/array.h
    include/alps/osiris/boost/ublas.h
    include/alps/osiris/buffer.h
    include/alps/osiris/comm.h
    include/alps/osiris/dump.h
    include/alps/osiris/dumparchive.h
    include/alps/osiris/mpdump.h
    include/alps/osiris/process.h
    include/alps/osiris/std/deque.h
    include/alps/osiris/std/impl.h
    include/alps/osiris/std/list.h
    include/alps/osiris/std/map.h
    include/alps/osiris/std/pair.h
    include/alps/osiris/std/set.h
    include/alps/osiris/std/stack.h
    include/alps/osiris/std/string.h
    include/alps/osiris/std/valarray.h
    include/alps/osiris/std/vector.h
    include/alps/osiris/xdrdump.h
    "${PROJECT_SOURCE_DIR}/third_party/xdr/include/alps/osiris/rpc_types.h"
    "${PROJECT_SOURCE_DIR}/third_party/xdr/include/alps/osiris/xdrcore.h")

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
