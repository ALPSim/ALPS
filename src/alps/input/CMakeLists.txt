# Application parameter parsers. Applications read their own parameters
# through these functions instead of from their task's Parameters directly.
target_sources(alps PRIVATE
  ctint.cpp
  dmft.cpp
  hirschfye.cpp
  hybridization.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET input_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/input/ctint.hpp
        include/alps/input/dmft.hpp
        include/alps/input/hirschfye.hpp
        include/alps/input/hybridization.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
