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

alps_register_module(NAME accumulators HEADER_SET accumulators_headers
  TARGETS alps
  DEPENDS alea containers hdf5 numerics utilities
  DEBT_FILE architecture-debt.json
  PRIVATE_FILES
    src/ngs/accumulator/accumulator.cpp
    src/ngs/accumulator/deprecated/alea/accumulator_set.cpp
    src/ngs/accumulator/deprecated/alea/result_set.cpp
)
