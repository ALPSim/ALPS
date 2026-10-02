# Utilities have no link dependency on the simulation or archive runtime.
add_library(alps_utilities
  src/copyright.cpp
  src/encode.cpp
  src/os.cpp
  src/short_print.cpp
  src/sleep.cpp
  src/stacktrace.cpp
  src/temporary_filename.cpp
  src/vmusage.cpp)
add_library(ALPS::utilities ALIAS alps_utilities)
set_target_properties(alps_utilities PROPERTIES
  EXPORT_NAME utilities
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN};$<TARGET_FILE_DIR:Boost::filesystem>")
target_link_libraries(alps_utilities PUBLIC alps_headers Boost::filesystem Threads::Threads)
alps_copy_runtime_dependencies(alps_utilities)
if(WIN32)
  target_link_libraries(alps_utilities PRIVATE ws2_32 mswsock)
endif()
generate_export_header(alps_utilities EXPORT_MACRO_NAME ALPS_UTILITIES_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/src/alps/utilities_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_UTILITIES_STATIC_DEFINE)
endif()
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  FILES "${PROJECT_BINARY_DIR}/src/alps/utilities_export.h")
target_sources(alps_headers INTERFACE
  FILE_SET utilities_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/cast.hpp
    include/alps/ngs/detail/remove_cvr.hpp
    include/alps/ngs/detail/type_wrapper.hpp
    include/alps/ngs/short_print.hpp
    include/alps/ngs/sleep.hpp
    include/alps/ngs/stacktrace.hpp
    include/alps/utility/assign.hpp
    include/alps/utility/bitops.hpp
    include/alps/utility/copyright.hpp
    include/alps/utility/data.hpp
    include/alps/utility/encode.hpp
    include/alps/utility/factory.hpp
    include/alps/utility/generator.hpp
    include/alps/utility/make_copy.hpp
    include/alps/utility/numeric_cast.hpp
    include/alps/utility/os.hpp
    include/alps/utility/resize.hpp
    include/alps/utility/set_zero.hpp
    include/alps/utility/size.hpp
    include/alps/utility/temporary_filename.hpp
    include/alps/utility/vectorio.hpp
    include/alps/utility/vmusage.hpp)

install(TARGETS alps_utilities EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
