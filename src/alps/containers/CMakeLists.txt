add_library(alps_containers INTERFACE)
add_library(ALPS::containers ALIAS alps_containers)
set_target_properties(alps_containers PROPERTIES EXPORT_NAME containers)
target_link_libraries(alps_containers INTERFACE ALPS::configuration)
# Containers source ownership.
target_sources(alps_containers INTERFACE
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
    include/alps/multi_array/io.hpp
    include/alps/multi_array/multi_array.hpp
    include/alps/multi_array/serialization.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
