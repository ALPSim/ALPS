# Compatibility with older consumers. New projects only need:
#   find_package(ALPS CONFIG REQUIRED)
#   target_link_libraries(my_target PRIVATE ALPS::alps)
include_guard(GLOBAL)
include_directories(${ALPS_INCLUDE_DIRS})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
include("${CMAKE_CURRENT_LIST_DIR}/add_alps_test.cmake")
# Never replace the consumer's compilers, flags or build configuration.
