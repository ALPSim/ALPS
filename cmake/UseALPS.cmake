# Compatibility entry point for projects using include(${ALPS_USE_FILE}).
# Usage requirements now come from the installed target instead of replacing
# the consumer's compiler and global flags with those of the SDK build.
include_guard(DIRECTORY)
if(NOT TARGET ALPS::alps)
  find_package(ALPS CONFIG REQUIRED)
endif()
link_libraries(ALPS::alps)
include("${CMAKE_CURRENT_LIST_DIR}/add_alps_test.cmake")
