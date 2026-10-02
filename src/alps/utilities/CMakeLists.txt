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
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_utilities PUBLIC ALPS::configuration ALPS::containers Boost::filesystem Threads::Threads)
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_utilities PUBLIC "$<COMPILE_ONLY:MPI::MPI_CXX;Boost::mpi>")
endif()
alps_dependency_rpath(alps_utilities)
alps_copy_runtime_dependencies(alps_utilities)
if(WIN32)
  target_link_libraries(alps_utilities PRIVATE ws2_32 mswsock)
endif()
generate_export_header(alps_utilities EXPORT_MACRO_NAME ALPS_UTILITIES_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/utilities_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_utilities PUBLIC ALPS_UTILITIES_STATIC_DEFINE)
endif()
target_sources(alps_utilities PUBLIC
  FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include"
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/utilities_export.h")
target_sources(alps_utilities PUBLIC
  FILE_SET utilities_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/cctype.h
    include/alps/factory.h
    include/alps/lambda.hpp
    include/alps/ngs/boost_mpi.hpp
    include/alps/ngs/cast.hpp
    include/alps/ngs/config.hpp
    include/alps/ngs/detail/remove_cvr.hpp
    include/alps/ngs/detail/type_wrapper.hpp
    include/alps/ngs/hash.hpp
    include/alps/ngs/mutex.hpp
    include/alps/ngs/short_print.hpp
    include/alps/ngs/sleep.hpp
    include/alps/ngs/stacktrace.hpp
    include/alps/ngs/stringify.hpp
    include/alps/ngs/thread_exceptions.hpp
    include/alps/stringvalue.h
    include/alps/type_traits/average_type.hpp
    include/alps/type_traits/change_value_type.hpp
    include/alps/type_traits/covariance_type.hpp
    include/alps/type_traits/element_type.hpp
    include/alps/type_traits/has_value_type.hpp
    include/alps/type_traits/is_complex.hpp
    include/alps/type_traits/is_scalar.hpp
    include/alps/type_traits/is_sequence.hpp
    include/alps/type_traits/is_symbolic.hpp
    include/alps/type_traits/iterator_type.hpp
    include/alps/type_traits/norm_type.hpp
    include/alps/type_traits/param_type.hpp
    include/alps/type_traits/real_type.hpp
    include/alps/type_traits/slice.hpp
    include/alps/type_traits/type_tag.hpp
    include/alps/utility/assign.hpp
    include/alps/utility/bitops.hpp
    include/alps/utility/copyright.hpp
    include/alps/utility/data.hpp
    include/alps/utility/encode.hpp
    include/alps/utility/factory.hpp
    include/alps/utility/generator.hpp
    include/alps/utility/make_copy.hpp
    include/alps/utility/os.hpp
    include/alps/utility/resize.hpp
    include/alps/utility/set_zero.hpp
    include/alps/utility/size.hpp
    include/alps/utility/temporary_filename.hpp
    include/alps/utility/vectorio.hpp
    include/alps/utility/vmusage.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME utilities HEADER_SET utilities_headers
  HEADER_TARGETS alps_utilities
  TARGETS alps alps_utilities
  DEPENDS configuration containers
  GENERATED_HEADERS alps/utilities_export.h
)
