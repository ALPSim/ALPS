# ALPSCore's parameters, alps::params in <alps/params.hpp>, copied from
# ALPSCore; see the README. They read INI files with the bundled iniparser C
# library, and have no link dependency on the legacy parameter/parser runtime.
add_library(alps_params_iniparser OBJECT iniparser/src/dictionary.c iniparser/src/iniparser.c)
set_target_properties(alps_params_iniparser PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(alps_params
  src/dict_value.cpp
  src/dictionary.cpp
  src/iniparser_interface.cpp
  src/params.cpp
  src/text_values.cpp
  $<TARGET_OBJECTS:alps_params_iniparser>)
target_include_directories(alps_params PRIVATE iniparser/src)
target_link_libraries(alps_params PUBLIC ALPS::configuration ALPS::utilities ALPS::hdf5)
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_params PUBLIC MPI::MPI_CXX)
endif()
alps_configure_component(alps_params params)
# ALPSCore's sources have no export macros; export everything as ALPSCore did.
set_target_properties(alps_params PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_sources(alps_params PUBLIC
  FILE_SET alpscore_params_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/dictionary.hpp
    include/alps/params.hpp
    include/alps/params/convenience_params.hpp
    include/alps/params/dict_exceptions.hpp
    include/alps/params/dict_types.hpp
    include/alps/params/dict_value.hpp
    include/alps/params/dict_value_impl.hpp
    include/alps/params/hdf5_variant.hpp
    include/alps/params/iniparser_interface.hpp
    include/alps/params/params_impl.hpp
    include/alps/params/serialize_variant.hpp
    # ALPS's addition
    include/alps/params/text_values.hpp)
if(ALPS_ENABLE_MPI)
  # It broadcasts through <alps/utilities/mpi_vector.hpp>, an MPI-only header.
  target_sources(alps_params PUBLIC FILE_SET alpscore_params_headers FILES
    include/alps/params/mpi_variant.hpp)
endif()
install(FILES COPYRIGHT.TXT LICENSE.TXT
  DESTINATION "${CMAKE_INSTALL_DATADIR}/alps/licenses/params" COMPONENT libraries)
install(FILES iniparser/LICENSE iniparser/AUTHORS
  DESTINATION "${CMAKE_INSTALL_DATADIR}/alps/licenses/params/iniparser" COMPONENT libraries)
add_subdirectory(adapters)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
