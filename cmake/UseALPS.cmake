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
    set(CMAKE_Fortran_COMPILER ${ALPS_CMAKE_Fortran_COMPILER} CACHE FILEPATH "Fortran compiler." FORCE)
    set(CMAKE_Fortran_FLAGS ${ALPS_CMAKE_Fortran_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG ${ALPS_CMAKE_Fortran_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE ${ALPS_CMAKE_Fortran_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
  endif(NOT PREVENT_ALPS_COMPILERS)

  # The Boost Root Dir used by ALPS
  set(Boost_ROOT_DIR ${ALPS_Boost_ROOT_DIR})
  set(Boost_INCLUDE_DIR ${ALPS_Boost_INCLUDE_DIR})
	set(Boost_LIBRARIES ${ALPS_Boost_LIBRARIES})

  # OpenMP
  set(OpenMP_C_FLAGS ${ALPS_OpenMP_C_FLAGS})
  set(OpenMP_CXX_FLAGS ${ALPS_OpenMP_CXX_FLAGS})

  # MPI
  set(MPI_FOUND ${ALPS_MPI_FOUND})
  set(MPI_DEFINITIONS ${ALPS_MPI_DEFINITONS})
  set(MPI_INCLUDE_DIR ${ALPS_MPI_INCLUDE_DIR})
  set(MPI_LINKER_FLAGS ${ALPS_MPI_LINKER_FLAGS})
  set(MPI_LIBRARIES ${ALPS_MPI_LIBRARIES})
  set(MPIEXEC ${ALPS_MPIEXEC})
  set(MPIEXEC_NUMPROC_FLAG ${ALPS_MPIEXEC_NUMPROC_FLAG})

  # BLAS/LAPACK
  set(LAPACK_FOUND ${ALPS_LAPACK_FOUND})
  set(LAPACK_DEFINITIONS ${ALPS_LAPACK_DEFINITIONS})
  set(LAPACK_LINKER_FLAGS ${ALPS_LAPACK_LINKER_FLAGS})
  set(LAPACK_LIBRARIES ${ALPS_LAPACK_LIBRARIES})
  set(LAPACK_LIBRARY ${ALPS_LAPACK_LIBRARY})
  set(BLAS_LIBRARIES ${ALPS_BLAS_LIBRARIES})
  set(BLAS_LIBRARY ${ALPS_BLAS_LIBRARY})
  set(MKL_INCLUDE_DIR ${ALPS_MKL_INCLUDE_DIR})

  # Python
  set(PYTHON_INTERPRETER ${ALPS_PYTHON_INTERPRETER})
  set(PYTHON_INCLUDE_DIRS ${ALPS_PYTHON_INCLUDE_DIRS})
  set(PYTHON_NUMPY_INCLUDE_DIR ${ALPS_PYTHON_NUMPY_INCLUDE_DIR})
  set(PYTHON_LIBRARY ${ALPS_PYTHON_LIBRARY})
  set(PYTHON_SITE_PKG ${ALPS_PYTHON_SITE_PKG})
  set(PYTHON_EXTRA_LIBS ${ALPS_PYTHON_EXTRA_LIBS})
  set(PYTHON_LINK_FOR_SHARED ${ALPS_PYTHON_LINK_FOR_SHARED})

  # FFTW
  set(FFTW_LIBRARIES ${ALPS_FFTW_LIBRARIES})
  set(FFTW_INCLUDE_DIR ${ALPS_FFTW_INCLUDE_DIR})

  # SQLite
  set(SQLite_FOUND ${ALPS_SQLite_FOUND})
  set(SQLite_INCLUDE_DIR ${ALPS_SQLite_INCLUDE_DIR})
  set(SQLite_LIBRARIES ${ALPS_SQLite_LIBRARIES})

  # LPSolve
  # set(LPSolve_FOUND ${ALPS_LPSolve_FOUND})
  # set(LPSolve_INCLUDE_DIR ${ALPS_LPSolve_INCLUDE_DIR})
  # set(LPSolve_LIBRARIES ${ALPS_LPSolve_LIBRARIES})

  # HDF5
  set(HDF5_FOUND ${ALPS_HDF5_FOUND})
  set(HDF5_INCLUDE_DIR ${ALPS_HDF5_INCLUDE_DIR})
  set(HDF5_LIBRARIES  ${ALPS_HDF5_LIBRARIES})
  set(HDF5_DIFF_EXECUTABLE ${ALPS_HDF5_DIFF_EXECUTABLE})
  
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

  # RPATH setting
  set(CMAKE_INSTALL_NAME_DIR "${ALPS_ROOT_DIR}/lib" FORCE)
  set(CMAKE_INSTALL_RPATH "${ALPS_ROOT_DIR}/lib" FORCE)

  # test macro
  include(${ALPS_ROOT_DIR}/share/alps/add_alps_test.cmake)
endif(NOT ALPS_USE_FILE_INCLUDED)
