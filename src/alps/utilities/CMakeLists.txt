# Utilities have no link dependency on the simulation or archive runtime.
add_library(alps_utilities
  src/copyright.cpp
  src/citations.cpp
  src/encode.cpp
  src/os.cpp
  src/short_print.cpp
  src/signal.cpp
  src/sleep.cpp
  src/stacktrace.cpp
  src/temporary_filename.cpp
  src/vmusage.cpp)
target_link_libraries(alps_utilities PUBLIC ALPS::configuration ALPS::containers Boost::filesystem Threads::Threads)

# Stacktraces require both platform headers and linkable runtime symbols.
# Detect them for every build, including embedded consumers on musl.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
block()
  cmake_push_check_state(RESET)
  set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
  check_cxx_source_compiles([[
    #include <cxxabi.h>
    #include <execinfo.h>
    #include <cstdlib>
    int main() {
      void *frames[1];
      int depth = backtrace(frames, 1);
      std::free(backtrace_symbols(frames, depth));
      int status = 0;
      std::free(abi::__cxa_demangle("", 0, 0, &status));
    }
  ]] ALPS_HAVE_STACKTRACE)
  cmake_pop_check_state()
endblock()
if(NOT ALPS_HAVE_STACKTRACE)
  target_compile_definitions(alps_utilities PRIVATE ALPS_NGS_NO_STACKTRACE)
endif()

if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_utilities PUBLIC "$<COMPILE_ONLY:MPI::MPI_CXX;Boost::mpi>")
endif()
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_utilities PUBLIC MPI::MPI_CXX Boost::mpi)
endif()
if(WIN32)
  target_link_libraries(alps_utilities PRIVATE ws2_32 mswsock)
endif()
alps_configure_component(alps_utilities utilities)
target_include_directories(alps_utilities PRIVATE "${PROJECT_BINARY_DIR}/generated/include")
target_sources(alps_utilities PUBLIC
  FILE_SET utilities_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/cctype.h
    include/alps/factory.h
    include/alps/lambda.hpp
    include/alps/ngs/boost_mpi.hpp
    include/alps/ngs/config.hpp
    include/alps/ngs/hash.hpp
    include/alps/ngs/mutex.hpp
    include/alps/ngs/sleep.hpp
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
    include/alps/utilities/cast.hpp
    include/alps/utilities/citations.hpp
    include/alps/utilities/copyright.hpp
    include/alps/utilities/data.hpp
    include/alps/utilities/encode.hpp
    include/alps/utilities/make_copy.hpp
    include/alps/utilities/os.hpp
    include/alps/utilities/remove_cvr.hpp
    include/alps/utilities/resize.hpp
    include/alps/utilities/short_print.hpp
    include/alps/utilities/signal.hpp
    include/alps/utilities/size.hpp
    include/alps/utilities/stacktrace.hpp
    include/alps/utilities/stringify.hpp
    include/alps/utilities/temporary_filename.hpp
    include/alps/utilities/type_wrapper.hpp
    include/alps/utilities/vectorio.hpp
    include/alps/utilities/vmusage.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
