# Containers source ownership.
target_sources(alps_headers INTERFACE
  FILE_SET containers_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/fixed_capacity/checking.h
    include/alps/fixed_capacity/deque_detail.h
    include/alps/fixed_capacity/uninitialized_array.h
    include/alps/fixed_capacity_deque.h
    include/alps/fixed_capacity_fwd.h
    include/alps/fixed_capacity_traits.h
    include/alps/fixed_capacity_vector.h
    include/alps/multi_array.hpp
    include/alps/multi_array/functions.hpp
    include/alps/multi_array/io.hpp
    include/alps/multi_array/multi_array.hpp
    include/alps/multi_array/operators.hpp
    include/alps/multi_array/serialization.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME containers HEADER_SET containers_headers
  TARGETS alps
  DEPENDS numerics
)
