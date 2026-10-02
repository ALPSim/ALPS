# Parapack source ownership.
target_sources(alps PRIVATE
  src/ngs/clone.cpp
  src/ngs/clone_info.cpp
  src/ngs/job.cpp
  src/ngs/parapack.cpp
  src/ngs/worker_factory.cpp
  src/parapack/clone.C
  src/parapack/clone_info.C
  src/parapack/filelock.C
  src/parapack/job.C
  src/parapack/logger.C
  src/parapack/mc_worker.C
  src/parapack/measurement.C
  src/parapack/option.C
  src/parapack/parapack.C
  src/parapack/process_impl.C
  src/parapack/queue.C
  src/parapack/rng_helper.C
  src/parapack/staging.C
  src/parapack/types.C
  src/parapack/util.C
  src/parapack/version.C
  src/parapack/worker_factory.C)

target_sources(alps_headers INTERFACE
  FILE_SET parapack_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/parapack/clone.h
    include/alps/ngs/parapack/clone_info.h
    include/alps/ngs/parapack/clone_info_p.h
    include/alps/ngs/parapack/clone_proxy.h
    include/alps/ngs/parapack/job.h
    include/alps/ngs/parapack/job_p.h
    include/alps/ngs/parapack/params_p.h
    include/alps/ngs/parapack/parapack.h
    include/alps/ngs/parapack/simulation_p.h
    include/alps/ngs/parapack/worker_factory.h
    include/alps/parapack/clone.h
    include/alps/parapack/clone_info.h
    include/alps/parapack/clone_info_p.h
    include/alps/parapack/clone_proxy.h
    include/alps/parapack/clone_timer.h
    include/alps/parapack/exchange.h
    include/alps/parapack/exchange_multi.h
    include/alps/parapack/exp_number.h
    include/alps/parapack/filelock.h
    include/alps/parapack/footprint.h
    include/alps/parapack/integer_range.h
    include/alps/parapack/job.h
    include/alps/parapack/job_p.h
    include/alps/parapack/logger.h
    include/alps/parapack/mc_worker.h
    include/alps/parapack/measurement.h
    include/alps/parapack/montecarlo.h
    include/alps/parapack/option.h
    include/alps/parapack/parapack.h
    include/alps/parapack/permutation.h
    include/alps/parapack/process.h
    include/alps/parapack/queue.h
    include/alps/parapack/rng_helper.h
    include/alps/parapack/simulation_p.h
    include/alps/parapack/staging.h
    include/alps/parapack/temperature_scan.h
    include/alps/parapack/types.h
    include/alps/parapack/util.h
    include/alps/parapack/version.h
    include/alps/parapack/wanglandau.h
    include/alps/parapack/worker.h
    include/alps/parapack/worker_factory.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME parapack HEADER_SET parapack_headers
  TARGETS alps
  DEPENDS alea configuration expression hdf5 lattice legacy_parameters mc model numerics osiris params random scheduler utilities xml
)
