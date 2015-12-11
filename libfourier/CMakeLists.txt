cmake_minimum_required(VERSION 3.0)

project (fourier)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
##we also need the blas classes
find_package(Greens REQUIRED)

include(${ALPS_USE_FILE})
include_directories(${GREENS_INCLUDE_DIR})

add_library(fourier
fourier_generator.cpp  fourier_transform.cpp)

target_link_libraries(fourier ${GREENS_LIBRARIES})
