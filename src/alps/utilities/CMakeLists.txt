# Module-owned implementation and public headers; linked into the common ALPS runtime.
target_sources(alps PRIVATE
  src/copyright.cpp
  src/encode.cpp
  src/os.cpp
  src/short_print.cpp
  src/sleep.cpp
  src/stacktrace.cpp
  src/temporary_filename.cpp
  src/vmusage.cpp)
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

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
