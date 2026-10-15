# Mc source ownership.
target_sources(alps PRIVATE
  src/ngs/api.cpp
  src/ngs/mcbase.cpp
  src/ngs/stop_callback.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET mc_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/check_schedule.hpp
    include/alps/mcbase.hpp
    include/alps/mcmpiadapter.hpp
    include/alps/ngs.hpp
    include/alps/ngs/api.hpp
    include/alps/progress_callback.hpp
    include/alps/stop_callback.hpp)

# Public Python integration spellings survive the move to Python ownership.
target_sources(alps_headers INTERFACE FILE_SET python_compatibility_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/detail/export_sim_to_python.hpp
        include/alps/python/save_observable_to_hdf5.hpp)

# ALPS::mc: the Monte Carlo module of the ALPS Core Libraries, copied from
# ALPSCore; see the README. It defines alps::mcbase, alps::mcmpiadapter,
# alps::check_schedule, alps::stop_callback and alps::random01, as the NGS
# sources above do, so for now it is built and tested in the tree only: it
# joins neither ALPS::alps nor the installed SDK.
add_library(alps_mc src/api.cpp src/mcbase.cpp src/stop_callback.cpp)
target_link_libraries(alps_mc PUBLIC ALPS::configuration ALPS::utilities ALPS::hdf5
  ALPS::params ALPS::accumulators)
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_mc PUBLIC MPI::MPI_CXX)
endif()
alps_configure_component(alps_mc mc)
# ALPSCore's sources have no export macros; export everything as ALPSCore did.
set_target_properties(alps_mc PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_sources(alps_mc PUBLIC
  FILE_SET alpscore_mc_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/mc/api.hpp
    include/alps/mc/check_schedule.hpp
    include/alps/mc/mcbase.hpp
    include/alps/mc/mpiadapter.hpp
    include/alps/mc/random01.hpp
    include/alps/mc/stop_callback.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
