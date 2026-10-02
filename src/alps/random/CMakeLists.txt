# Random source ownership.
target_sources(alps PRIVATE
  src/random/parallel/get_prime.cpp
  src/random/parallel/primelist_64.cpp
  src/random/rngfactory.C)

target_sources(alps_headers INTERFACE
  FILE_SET random_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/random01.hpp
    include/alps/random.h
    include/alps/random/buffered_generator.hpp
    include/alps/random/buffered_rng.h
    include/alps/random/buffered_uniform_01.hpp
    include/alps/random/mersenne_twister.hpp
    include/alps/random/multivariate_normal_distribution.hpp
    include/alps/random/parallel.hpp
    include/alps/random/parallel/detail/get_prime.hpp
    include/alps/random/parallel/detail/primelist_64.hpp
    include/alps/random/parallel/detail/seed_macros.hpp
    include/alps/random/parallel/keyword.hpp
    include/alps/random/parallel/lcg64.hpp
    include/alps/random/parallel/mersenne_twister.hpp
    include/alps/random/parallel/mpi.hpp
    include/alps/random/parallel/primelist_64.hpp
    include/alps/random/parallel/seed.hpp
    include/alps/random/parallel/well.hpp
    include/alps/random/pseudo_des.h
    include/alps/random/random_choice.hpp
    include/alps/random/rngfactory.h
    include/alps/random/seed.h
    include/alps/random/sprng.hpp
    include/alps/random/sprng/cmrg.hpp
    include/alps/random/sprng/detail/buffer.hpp
    include/alps/random/sprng/detail/implementation.hpp
    include/alps/random/sprng/detail/interface.hpp
    include/alps/random/sprng/keyword.hpp
    include/alps/random/sprng/lcg.hpp
    include/alps/random/sprng/lcg64.hpp
    include/alps/random/sprng/lfg.hpp
    include/alps/random/sprng/mlfg.hpp
    include/alps/random/sprng/pmlcg.hpp
    include/alps/random/uniform_on_sphere_n.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME random HEADER_SET random_headers
  TARGETS alps
  DEPENDS configuration hdf5 utilities
)
