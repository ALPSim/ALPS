# ALPS::mc: the Monte Carlo module of the ALPS Core Libraries, copied from
# ALPSCore; see the README. It is part of ALPS::alps.
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
    include/alps/mc/stop_callback.hpp
    # ALPS's addition
    include/alps/progress_callback.hpp)
install(FILES ACKNOWLEDGE.TXT COPYRIGHT.TXT LICENSE.TXT
  DESTINATION "${CMAKE_INSTALL_DATADIR}/alps/licenses/mc" COMPONENT libraries)

# Public Python integration spellings survive the move to Python ownership.
target_sources(alps_headers INTERFACE FILE_SET python_compatibility_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/detail/export_sim_to_python.hpp
        include/alps/python/save_observable_to_hdf5.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
