# Application parameter parsers. Applications read their own parameters
# through these functions instead of from their task's Parameters directly.
target_sources(alps PRIVATE
  ctint.cpp
  hybridization.cpp
  looper.cpp
  parapack.cpp
  simplemc.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET input_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/input/ctint.hpp
        include/alps/input/hybridization.hpp
        include/alps/input/looper.hpp
        include/alps/input/parapack.hpp
        include/alps/input/simplemc.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
