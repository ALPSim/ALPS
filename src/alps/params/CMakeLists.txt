# Typed parameters have no link dependency on the legacy parameter/parser runtime.
add_library(alps_params
  src/ngs/paramproxy.cpp
  src/ngs/params.cpp
  src/ngs/paramvalue.cpp)
target_link_libraries(alps_params PUBLIC ALPS::configuration ALPS::hdf5 Boost::serialization
  PRIVATE "$<COMPILE_ONLY:ALPS::numerics>")
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_params PUBLIC MPI::MPI_CXX Boost::mpi)
endif()
alps_configure_component(alps_params params)
target_sources(alps_params PUBLIC
  FILE_SET params_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/detail/paramiterator.hpp
    include/alps/ngs/detail/paramproxy.hpp
    include/alps/ngs/detail/params_impl_base.hpp
    include/alps/ngs/detail/paramvalue.hpp
    include/alps/ngs/detail/paramvalue_reader.hpp
    include/alps/ngs/params.hpp)

add_subdirectory(adapters)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
