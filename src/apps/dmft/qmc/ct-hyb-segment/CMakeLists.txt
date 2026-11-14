cmake_minimum_required(VERSION 3.12)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# Not strictly necessary, but will create executables that are
# aware of ALPSCore location
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

project(alps_cthyb CXX)

# Disable in-source builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Request the ALPSCore package:
find_package(ALPSCore REQUIRED COMPONENTS accumulators params mc hdf5 gf)
find_package(Boost REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI)
find_package(NFFT3 REQUIRED)

set (ALPS_CTHYB_SOURCE
   hybsim.cpp
   hyblocal.cpp
   hybint.cpp
   hybfun.cpp
   hybretintfun.cpp
   hybmatrix.cpp
   hybmatrix_nfft.cpp
   hybconfig.cpp
   hybupdates.cpp
   hybevaluate.cpp
   hybmeasurements.cpp)
add_library(alps_hyb_exp ${ALPS_CTHYB_SOURCE})
target_compile_definitions(alps_hyb_exp
    PUBLIC
        BOOST_DISABLE_ASSERTS
        $<$<NOT:$<CONFIG:Debug>>:NDEBUG>
)
target_include_directories(alps_hyb_exp
    PUBLIC
        "${CMAKE_SOURCE_DIR}"
        "${CMAKE_BINARY_DIR}"
        ${Boost_INCLUDE_DIRS}
        ${MPI_CXX_INCLUDE_DIRS}
        ${NFFT3_INCLUDE_DIRS}
)
target_link_libraries(alps_hyb_exp
    PUBLIC
        ${ALPSCore_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        ${NFFT3_LIBRARIES}
)
set_target_properties(alps_hyb_exp PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(EXEC_NAME alps_cthyb)
add_executable(${EXEC_NAME} hybmain.cpp)
target_link_libraries(${EXEC_NAME}
    PRIVATE
        alps_hyb_exp
        ${ALPSCore_LIBRARIES}
        ${LAPACK_LIBRARIES}
)
set_target_properties(${EXEC_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
install(TARGETS ${EXEC_NAME} DESTINATION bin)

### Additional installation files
# configuration file for module managers, such as LMOD
configure_file("${CMAKE_SOURCE_DIR}/misc/alps_cthyb.lmod.in" "${CMAKE_BINARY_DIR}/alps_cthyb.lmod")
