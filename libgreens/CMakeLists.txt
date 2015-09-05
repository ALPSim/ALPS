cmake_minimum_required(VERSION 2.8)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

#we also need the blas classes
include_directories($ENV{BLAS_CLASSES_DIR})

add_library (green green_function.cpp nambu_green_function.cpp)
