# Accumulators source ownership.
if(ALPS_NGS_USE_NEW_ALEA)
  target_sources(alps PRIVATE src/ngs/accumulator/accumulator.cpp)
endif()

target_sources(alps_headers INTERFACE
  FILE_SET accumulators_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/accumulator.hpp
    include/alps/ngs/accumulator/accumulator.hpp
    include/alps/ngs/accumulator/deprecated/accumulator.hpp
    include/alps/ngs/accumulator/deprecated/accumulator/accumulator_impl.hpp
    include/alps/ngs/accumulator/deprecated/accumulator/arguments.hpp
    include/alps/ngs/accumulator/deprecated/accumulator_set.hpp
    include/alps/ngs/accumulator/deprecated/alea.hpp
    include/alps/ngs/accumulator/deprecated/feature/autocorrelation.hpp
    include/alps/ngs/accumulator/deprecated/feature/converged.hpp
    include/alps/ngs/accumulator/deprecated/feature/error.hpp
    include/alps/ngs/accumulator/deprecated/feature/feature_traits.hpp
    include/alps/ngs/accumulator/deprecated/feature/features.hpp
    include/alps/ngs/accumulator/deprecated/feature/fixed_size_binning.hpp
    include/alps/ngs/accumulator/deprecated/feature/generate_property.hpp
    include/alps/ngs/accumulator/deprecated/feature/histogram.hpp
    include/alps/ngs/accumulator/deprecated/feature/log_binning.hpp
    include/alps/ngs/accumulator/deprecated/feature/max_num_binning.hpp
    include/alps/ngs/accumulator/deprecated/feature/mean.hpp
    include/alps/ngs/accumulator/deprecated/feature/tags.hpp
    include/alps/ngs/accumulator/deprecated/feature/tau.hpp
    include/alps/ngs/accumulator/deprecated/feature/value_type.hpp
    include/alps/ngs/accumulator/deprecated/feature/weight.hpp
    include/alps/ngs/accumulator/deprecated/features.hpp
    include/alps/ngs/accumulator/deprecated/result.hpp
    include/alps/ngs/accumulator/deprecated/result_set.hpp
    include/alps/ngs/accumulator/deprecated/wrapper/accumulator_wrapper.hpp
    include/alps/ngs/accumulator/deprecated/wrapper/base_wrapper.hpp
    include/alps/ngs/accumulator/deprecated/wrapper/derived_wrapper.hpp
    include/alps/ngs/accumulator/deprecated/wrapper/result_type_wrapper.hpp
    include/alps/ngs/accumulator/deprecated/wrapper/result_wrapper.hpp
    include/alps/ngs/accumulator/feature.hpp
    include/alps/ngs/accumulator/feature/binning_analysis.hpp
    include/alps/ngs/accumulator/feature/count.hpp
    include/alps/ngs/accumulator/feature/error.hpp
    include/alps/ngs/accumulator/feature/max_num_binning.hpp
    include/alps/ngs/accumulator/feature/mean.hpp
    include/alps/ngs/accumulator/feature/weight.hpp
    include/alps/ngs/accumulator/feature/weight_holder.hpp
    include/alps/ngs/accumulator/parameter.hpp
    include/alps/ngs/accumulator/wrappers.hpp
    include/alps/ngs/mpi.hpp
    include/alps/ngs/result.hpp)

# ALPS::accumulators: the accumulators library of the ALPS Core Libraries,
# copied from ALPSCore. README.md records the source commit, the local changes
# and the license. It is separate from the NGS accumulators above, which use
# <alps/ngs/accumulator/...> and namespace alps::accumulator.
add_library(alps_accumulators
  src/accumulator_wrapper.cpp
  src/mpi.cpp
  src/namedaccumulators.cpp
  src/result_wrapper.cpp
  src/wrapper_set.cpp
  src/wrapper_set_hdf5.cpp
  src/feature/binning_analysis.cpp
  src/feature/count.cpp
  src/feature/error.cpp
  src/feature/max_num_binning.cpp
  src/feature/mean.cpp)
# Only numerics' headers: the accumulators need its vector arithmetic, not BLAS.
target_link_libraries(alps_accumulators PUBLIC ALPS::configuration ALPS::utilities ALPS::hdf5
  "$<COMPILE_ONLY:ALPS::numerics>")
if(ALPS_ENABLE_MPI)
  # It reduces through ALPSCore's alps/utilities/mpi.hpp and boost_mpi.hpp.
  target_link_libraries(alps_accumulators PUBLIC MPI::MPI_CXX)
endif()
alps_configure_component(alps_accumulators accumulators)
# ALPSCore's accumulators have no export macros; export everything as ALPSCore did.
set_target_properties(alps_accumulators PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_sources(alps_accumulators PUBLIC
  FILE_SET alpscore_accumulators_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/accumulators.hpp
    include/alps/accumulators/accumulator.hpp
    include/alps/accumulators/archive_traits.hpp
    include/alps/accumulators/convergence.hpp
    include/alps/accumulators/feature.hpp
    include/alps/accumulators/feature/binning_analysis.hpp
    include/alps/accumulators/feature/count.hpp
    include/alps/accumulators/feature/error.hpp
    include/alps/accumulators/feature/max_num_binning.hpp
    include/alps/accumulators/feature/mean.hpp
    include/alps/accumulators/feature/weight.hpp
    include/alps/accumulators/feature/weight_holder.hpp
    include/alps/accumulators/mpi.hpp
    include/alps/accumulators/namedaccumulators.hpp
    include/alps/accumulators/parameter.hpp
    include/alps/accumulators/wrapper_set.hpp
    include/alps/accumulators/wrappers.hpp)
install(FILES ACKNOWLEDGE.TXT COPYRIGHT.TXT LICENSE.TXT
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/accumulators
  COMPONENT libraries)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
