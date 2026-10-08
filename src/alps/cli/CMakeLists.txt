# Existing command-line grammars are independent of simulation and input parsing.
add_library(alps_cli src/mcoptions.cpp src/parseargs.cpp src/cli.cpp)
target_link_libraries(alps_cli PUBLIC ALPS::configuration ALPS::utilities Boost::program_options)
alps_configure_component(alps_cli cli)
target_sources(alps_cli PUBLIC
  FILE_SET cli_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/mcoptions.hpp include/alps/parseargs.hpp
        include/alps/utility/cli.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
