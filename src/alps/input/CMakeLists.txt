# Application parameter parsers. Applications read their own parameters
# through these functions instead of from their task's Parameters directly.
target_sources(alps PRIVATE
  dmft.cpp
  hirschfye.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET input_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/input/dmft.hpp
        include/alps/input/hirschfye.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
