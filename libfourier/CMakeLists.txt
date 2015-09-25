cmake_minimum_required(VERSION 3.0)

project (fourier)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need ALPSCore and the GFTools package
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package(GFTools REQUIRED)

include_directories(${GFTOOLS_INCLUDE_DIR})

add_library(fourier
fourier_generator.cpp  fourier_transform.cpp)

target_link_libraries(fourier ${GFTOOLS_LIBRARIES} ${ALPSCore_LIBRARIES})
