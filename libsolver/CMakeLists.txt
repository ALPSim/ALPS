cmake_minimum_required(VERSION 3.1)

project (solver)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

#find_package(ALPSCore REQUIRED COMPONENTS hdf5 params)


#add_library(solver external_solver.cpp)

#target_link_libraries(solver ${GFTOOLS_LIBRARIES} ${ALPSCore_LIBRARIES})
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} DESTINATION include)