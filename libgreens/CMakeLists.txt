cmake_minimum_required(VERSION 3.0)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
#we also need the blas classes
find_package(BlasClasses REQUIRED)

include(${ALPS_USE_FILE})
include_directories(${BLAS_CLASSES_INCLUDE_DIR})

add_library (green green_function.cpp nambu_green_function.cpp)
