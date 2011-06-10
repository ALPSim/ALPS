#
# This module is provided as ALPS_USE_FILE by ALPSConfig.cmake.  It can
# be INCLUDEd in a project to load the needed compiler and linker
# settings to use ALPS.
#

if(NOT ALPS_USE_FILE_INCLUDED)
  set(ALPS_USE_FILE_INCLUDED 1)

  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${ALPS_BUILD_TYPE} CACHE STRING "Type of build" FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
  set(BUILD_SHARED_LIBS ${ALPS_BUILD_SHARED_LIBS})
  list(APPEND CMAKE_MODULE_PATH ${ALPS_ROOT_DIR}/share/alps)

  # compilers and common options
  if(NOT PREVENT_ALPS_COMPILERS)
    set(CMAKE_C_COMPILER ${ALPS_CMAKE_C_COMPILER} CACHE FILEPATH "C compiler." FORCE)
    set(CMAKE_C_FLAGS ${ALPS_CMAKE_C_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS_DEBUG ${ALPS_CMAKE_C_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_C_FLAGS_RELEASE ${ALPS_CMAKE_C_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_CXX_COMPILER ${ALPS_CMAKE_CXX_COMPILER} CACHE FILEPATH "CXX compiler." FORCE)
    set(CMAKE_CXX_FLAGS ${ALPS_CMAKE_CXX_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_CXX_FLAGS_DEBUG ${ALPS_CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE ${ALPS_CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
  endif(NOT PREVENT_ALPS_COMPILERS)

  # Add macro definitions needed to use ALPS and dependent libraries
  add_definitions(${ALPS_EXTRA_DEFINITIONS})

  # Add include directories needed to use ALPS and dependent libraries
  include_directories(${ALPS_INCLUDE_DIRS} ${ALPS_EXTRA_INCLUDE_DIRS})

  # Add link directories needed to use ALPS and dependent libraries
  link_directories(${ALPS_LIBRARY_DIRS})

  # Add linker flags needed to use ALPS and dependent libraries
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ALPS_EXTRA_LINKER_FLAGS} ${ALPS_CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "Flags used by the linker." FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${ALPS_CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING "Flags used by the linker during debug builds." FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${ALPS_CMAKE_EXE_LINKER_FLAGS_RELEASE}" CACHE STRING "Flags used by the linker during release builds." FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ALPS_EXTRA_LINKER_FLAGS} ${ALPS_CMAKE_EXE_CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING "Flags used by the linker during the creation of dll's." FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${ALPS_CMAKE_SHARED_LINKER_FLAGS_DEBUG}" CACHE STRING "Flags used by the linker during debug builds." FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} ${ALPS_CMAKE_SHARED_LINKER_FLAGS_RELEASE}" CACHE STRING "Flags used by the linker during release builds." FORCE)

  # list of ALPS and dependent libraries
  set(ALPS_LIBRARIES alps boost ${ALPS_EXTRA_LIBRARIES} CACHE STRING "List of ALPS and dependent libraries." FORCE)

  # RPATH setting
  set(CMAKE_INSTALL_NAME_DIR "${ALPS_ROOT_DIR}/lib" FORCE)
  set(CMAKE_INSTALL_RPATH "${ALPS_ROOT_DIR}/lib" FORCE)

  # test macro
  include(${ALPS_ROOT_DIR}/share/alps/add_alps_test.cmake)
endif(NOT ALPS_USE_FILE_INCLUDED)
