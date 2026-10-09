# Accumulators source ownership.

target_sources(alps_headers INTERFACE
  FILE_SET accumulators_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/accumulator.hpp
    include/alps/ngs/accumulator/accumulator.hpp
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
