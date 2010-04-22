#
# This module is provided as ALPS_USE_FILE by ALPSConfig.cmake.  It can
# be INCLUDEd in a project to load the needed compiler and linker
# settings to use ALPS.
#

if(NOT ALPS_USE_FILE_INCLUDED)
  set(ALPS_USE_FILE_INCLUDED 1)

  set(CMAKE_BUILD_TYPE ${ALPS_BUILD_TYPE} CACHE STRING "Type of build" FORCE)
  set(BUILD_SHARED_LIBS ${ALPS_BUILD_SHARED_LIBS})

  # compilers and common options
  set(CMAKE_C_COMPILER ${ALPS_CMAKE_C_COMPILER})
  set(CMAKE_C_FLAGS ${ALPS_CMAKE_C_FLAGS})
  set(CMAKE_C_FLAGS_DEBUG ${ALPS_CMAKE_C_FLAGS_DEBUG})
  set(CMAKE_C_FLAGS_RELEASE ${ALPS_CMAKE_C_FLAGS_RELEASE})
  set(CMAKE_CXX_COMPILER ${ALPS_CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS ${ALPS_CMAKE_CXX_FLAGS})
  set(CMAKE_CXX_FLAGS_DEBUG ${ALPS_CMAKE_CXX_FLAGS_DEBUG})
  set(CMAKE_CXX_FLAGS_RELEASE ${ALPS_CMAKE_CXX_FLAGS_RELEASE})

  # Add macro definitions needed to use ALPS and dependent libraries
  add_definitions(${ALPS_EXTRA_DEFINITIONS})

  # Add include directories needed to use ALPS and dependent libraries
  include_directories(${ALPS_INCLUDE_DIRS} ${ALPS_EXTRA_INCLUDE_DIRS})

  # Add link directories needed to use ALPS and dependent libraries
  link_directories(${ALPS_LIBRARY_DIRS})

  # Add linker flags needed to use ALPS and dependent libraries
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ALPS_EXTRA_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ALPS_EXTRA_LINKER_FLAGS}")

  # list of ALPS and dependent libraries
  set(ALPS_LIBRARIES alps boost ${ALPS_EXTRA_LIBRARIES})

  # test macro
  include(${CMAKE_INSTALL_PREFIX}/share/alps/add_alps_test.cmake)
endif(NOT ALPS_USE_FILE_INCLUDED)
