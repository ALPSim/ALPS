# - FindLapack.cmake
# Find LAPACK/BLAS (and compatible) numerical libraries
#
# This module will define the following values:
# LAPACK_LIBRARY
# BLAS_LIBRARY
#
# 0) check if BLAS_LIBRARY and LAPACK_LIBRARY have already been set on command line, if so: take these values.
# 1) use ENV MKL
# 2) if Intel compiler pick MKL (-mkl)
# 3) search MKL in usual paths
# 4) search ENV ATLAS
# 5) search generic lapack/blas using CMake-provided FindLAPACK.cmake
#     if BLA_VENDOR is set, try to use that vendor's LAPACK implementation
# 6) if build is on cray use hardcoded path
# 7) give up

#  Copyright (C)  2009-2010 Matthias Troyer <troyer@comp-phys.org>
#  Copyright (C)  2009-2010 Synge Todo <wistaria@comp-phys.org>
#  Copyright (C)  2009-2010 Bela Bauer
#  Copyright (C)  2009-2010 Brigitte Surer
#  Copyright (C)  2009-2010 Lukas Gamper
#  Copyright (C)  2009-2012 Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
#  Copyright (C)       2010 Emanuel Gull <gull@phys.columbia.edu>
#  Copyright (C)       2012 Michele Dolfi <dolfim@phys.ethz.ch>
#
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# TODO: - look for MKL on Windows
#       - look for goto2

# Includes
include(CheckFunctionExists)


###################################################
# BLAS_LIBRARY, LAPACK_LIBRARY manually defined.
###################################################

IF(BLAS_LIBRARY AND LAPACK_LIBRARY)
 SET(LAPACK_LIBRARY_INIT 1)
 SET(BLAS_LIBRARY_INIT 1)
ENDIF(BLAS_LIBRARY AND LAPACK_LIBRARY)

########################################################################################################
# Looking for MKL.
#   For parallel MKL, OpenMP check has to be done beforehand.
# 0) $ENV{MKL} can be defined from http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#    if specified, this settings are chosen
# 1) If compiler is Intel >= 12 (Intel Composer XE 2011/2013, Intel Compiler Pro)
# 1.1) When OPENMP_FOUND=ON and ALPS_USE_MKL_PARALLEL=ON, use -mkl=parallel
# 1.2) When OPENMP_FOUND=OFF or ALPS_USE_MKL_PARALLEL=OFF, use -mkl=sequential
# 2) If $ENV{MKLROOT} / $ENV{MKL_HOME} defined (done by MKL tools/environment scripts), use the linking from advisor
# 3) Look for MKL libraries in MKL_PATHS
########################################################################################################

SET(MKL_PATHS "/usr/local/lib /usr/lib")

# 0) $ENV{MKL} can be set for explicit linking
if($ENV{MKL} MATCHES "mkl")
  set(LAPACK_LIBRARY $ENV{MKL})
  set(BLAS_LIBRARY "")
  set(HAVE_MKL TRUE)
endif($ENV{MKL} MATCHES "mkl")

# 1) Intel compiler >= 12
if(NOT HAVE_MKL AND NOT LAPACK_LIBRARY_INIT)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(INTEL_WITH_MKL FALSE)
    if(DEFINED CMAKE_CXX_COMPILER_VERSION)
      if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 12
         OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 12)
        set(INTEL_WITH_MKL TRUE)
      endif()
    endif()
    if($ENV{MKLROOT} MATCHES "composer")
      set(INTEL_WITH_MKL TRUE)
    endif()
    
    if(INTEL_WITH_MKL)
      if(OPENMP_FOUND AND ALPS_USE_MKL_PARALLEL)
        set(LAPACK_LIBRARY "-mkl=parallel")
      else(OPENMP_FOUND AND ALPS_USE_MKL_PARALLEL)
        set(LAPACK_LIBRARY "-mkl=sequential")
      endif(OPENMP_FOUND AND ALPS_USE_MKL_PARALLEL)
      set(BLAS_LIBRARY "")
      set(HAVE_MKL TRUE)
    endif(INTEL_WITH_MKL)
  endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
endif(NOT HAVE_MKL AND NOT LAPACK_LIBRARY_INIT)

# 2) Use MKLROOT / MKL_HOME and standard linking
if(NOT HAVE_MKL AND NOT LAPACK_LIBRARY_INIT)
  set(mkl_home "")
  if($ENV{MKLROOT} MATCHES "mkl")
    set(mkl_home $ENV{MKLROOT})
  elseif($ENV{MKL_HOME} MATCHES "mkl")
    set(mkl_home $ENV{MKL_HOME})
  endif()

  if(mkl_home MATCHES "mkl")
    file( STRINGS "${mkl_home}/include/mkl.h" _mkl_h_content REGEX "__INTEL_MKL" )
    string(REGEX REPLACE ".*#define __INTEL_MKL__ ([0-9]+).*"        "\\1" MKL_VERSION_MAJOR  "${_mkl_h_content}")
    string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__ ([0-9]+).*"  "\\1" MKL_VERSION_MINOR  "${_mkl_h_content}")
    string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__ ([0-9]+).*" "\\1" MKL_VERSION_UPDATE "${_mkl_h_content}")
    set(MKL_VERSION "${MKL_VERSION_MAJOR}.${MKL_VERSION_MINOR}.${MKL_VERSION_UPDATE}")
    
    # STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" MKL_VERSION ${mkl_home})
    # set(MKL_VERSION_RAW ${MKL_VERSION})
    # if((DEFINED MKL_VERSION OR MKL_VERSION VERSION_LESS 10 OR MKL_VERSION VERSION_GREATER 2000) AND ${mkl_home} MATCHES "composer")
    #   set(MKL_VERSION "10.3.1") # set to arbitrary 10.3 version (they should behave all the same)
    # endif((DEFINED MKL_VERSION OR MKL_VERSION VERSION_LESS 10 OR MKL_VERSION VERSION_GREATER 2000) AND ${mkl_home} MATCHES "composer")
    # if((DEFINED MKL_VERSION OR MKL_VERSION VERSION_LESS 10 OR MKL_VERSION VERSION_GREATER 2000) AND ${mkl_home} MATCHES "Compiler")
    #   set(MKL_VERSION "10.2.4") # set to arbitrary 10.2 version (they should behave all the same)
    # endif((DEFINED MKL_VERSION OR MKL_VERSION VERSION_LESS 10 OR MKL_VERSION VERSION_GREATER 2000) AND ${mkl_home} MATCHES "Compiler")
    
    message(STATUS "LAPACK DEBUG::Compiler id ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "LAPACK DEBUG::Compiler version ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "LAPACK DEBUG::Processor ${CMAKE_SYSTEM_PROCESSOR}")
    message(STATUS "LAPACK DEBUG::ENV{MKL} $ENV{MKL}")
    message(STATUS "LAPACK DEBUG::ENV{MKLROOT} $ENV{MKLROOT}")
    message(STATUS "LAPACK DEBUG::ENV{MKL_HOME} $ENV{MKL_HOME}")
    # message(STATUS "LAPACK DEBUG::MKL_VERSION (raw) ${MKL_VERSION_RAW}")
    message(STATUS "LAPACK DEBUG::MKL_VERSION ${MKL_VERSION}")
    
    # OS thread library (pthread required by MKL)
    find_package(Threads REQUIRED)
    # MKL core
    if(OPENMP_FOUND AND ALPS_USE_MKL_PARALLEL)
      # No parallel mode support for MKL < 10.0
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # Intel with Intel OpenMP
        set(MKL_CORE -lmkl_intel_thread -lmkl_core -liomp5)
      elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        # GCC with GNU OpenMP
        # MKL with g++ needs gfortran
        set(MKL_CORE -lmkl_gnu_thread -lmkl_core -lgfortran)
      endif()
    else()
      if(${MKL_VERSION} MATCHES "1[0-1]\\.[0-3]\\.[0-9]+")
        set(MKL_CORE -lmkl_sequential -lmkl_core)
      else() # MKL < 10.0
        set(MKL_CORE -lmkl_lapack -lmkl -lguide)
      endif()
    endif()
    # basic data type model interface
    # - assuming ILP32 or LP64
    if(${MKL_VERSION} MATCHES "1[0-1]\\.[0-3]\\.[0-9]+")
      if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
        set(MKL_INTERFACE -lmkl_intel_lp64)
      elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686")
        set(MKL_INTERFACE -lmkl_intel)
      else()
        message(SEND_ERROR "MKL: the processor type of this system is not supported")
      endif()
    else() # MKL < 10.0
      set(MKL_INTERFACE "")
    endif()
    # MKL library path
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      if(${MKL_VERSION} MATCHES "11\\.[0-9]\\.[0-9]+" OR ${MKL_VERSION} MATCHES "10\\.3\\.[0-9]+")
        set(MKL_LIBRARY_PATH -L${mkl_home}/lib)
      else() # MKL < 10.3
        if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/em64t)
        elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/32)
        else()
          message(SEND_ERROR "MKL: the processor type of this system is not supported")
        endif()
      endif()
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      if(${MKL_VERSION} MATCHES "11\\.[0-9]\\.[0-9]+" OR ${MKL_VERSION} MATCHES "10\\.3\\.[0-9]+")
        if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/intel64)
        elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/ia32)
        else()
          message(SEND_ERROR "MKL: the processor type of this system is not supported")
        endif()
      else() # MKL < 10.3 have the same PATH
        if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/em64t)
        elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/32)
        elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
          set(MKL_LIBRARY_PATH -L${mkl_home}/lib/64)
        else()
          message(SEND_ERROR "MKL: the processor type of this system is not supported")
        endif()
      endif()
    endif()
    # combine together
    set(LAPACK_LIBRARY ${MKL_LIBRARY_PATH} ${MKL_INTERFACE} ${MKL_CORE} ${CMAKE_THREAD_LIBS_INIT} -lm)

    # unset local variables
    unset(MKL_LIBRARY_PATH)
    unset(MKL_INTERFACE)
    unset(MKL_CORE)

    set(BLAS_LIBRARY "")
    set(HAVE_MKL TRUE)
  endif(mkl_home MATCHES "mkl")
endif(NOT HAVE_MKL AND NOT LAPACK_LIBRARY_INIT)

IF(HAVE_MKL)
  # Checking if it works
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARY} ${LAPACK_LIBRARY})
  check_function_exists("sgemm" _libraries_work)
  set(CMAKE_REQUIRED_LIBRARIES)
  
  if(NOT _libraries_work)
    message(WARNING "MKL was detected but I'm not able to use it.")
    message(STATUS "MKL settings were:")
    message(STATUS "   BLAS_LIBRARY = ${BLAS_LIBRARY}")
    message(STATUS "   LAPACK_LIBRARY = ${LAPACK_LIBRARY}")
    set(BLAS_LIBRARY)
    set(LAPACK_LIBRARY)
    set(HAVE_MKL)
  else()
    message(STATUS "Found intel/mkl library")
    set(LAPACK_LIBRARY_INIT 1)
    set(BLAS_LIBRARY_INIT 1)
    set(MKL_INC_PATHS $ENV{mkl_home}/include ${MKL_PATHS}) 
    find_path(MKL_INCLUDE_DIR mkl.h ${MKL_INC_PATHS})
    include_directories(${MKL_INCLUDE_DIR})
    set(ALPS_HAVE_MKL 1) # MKL flag set in alps/config.h
  endif(NOT _libraries_work)
ENDIF(HAVE_MKL)


# SET(MKL_PATHS 
#   $ENV{MKL_HOME}/lib/intel64
#   $ENV{MKL_HOME}/lib/em${QMC_BITS}t
#   $ENV{MKL_HOME}/lib/${QMC_BITS}
#   $ENV{MKL_HOME}/lib
#   ${MKL_PATHS} 
#   ) 
# #MESSAGE(STATUS "Looking for intel/mkl library in ${MKL_PATHS}")
# 
# IF(NOT LAPACK_LIBRARY_INIT)
#   IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
#     IF(LINK_LAPACK_OLD)
#       FIND_LIBRARY(LAPACK_LIBRARY NAMES mkl_lapack PATHS ${MKL_PATHS})
#       FIND_LIBRARY(BLAS_LIBRARY NAMES mkl PATHS ${MKL_PATHS})
#     ELSE(LINK_LAPACK_OLD)
#       FIND_LIBRARY(LAPACK_LIBRARY STATIC NAMES mkl_lapack PATHS ${MKL_PATHS})
#       FIND_LIBRARY(BLAS_LIBRARY  NAMES mkl_em64t PATHS ${MKL_PATHS})
#       MESSAGE("-- mkl 10.0.[0-2] warning for EM64T")
#       MESSAGE("-- Replace libmkl_lapack.so in CMakeCache.txt by libmkl_lapack.a")
#     ENDIF(LINK_LAPACK_OLD)
#   ELSE(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
#     FIND_LIBRARY(LAPACK_LIBRARY 
#       NAMES mkl_lapack 
#       PATHS ${MKL_PATHS}
#       )
#     FIND_LIBRARY(BLAS_LIBRARY
#       NAMES mkl
#       PATHS ${MKL_PATHS}
#     )
#   ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
# 
#   FIND_LIBRARY(INTEL_GUIDE_LIBRARY
#     NAMES guide
#     PATHS ${MKL_PATHS}
#   )
#   MARK_AS_ADVANCED(INTEL_GUIDE_LIBRARY)
#   IF(NOT INTEL_GUIDE_LIBRARY MATCHES "NOTFOUND")
#     SET(REQUIRE_PTHREAD TRUE)
#     FIND_LIBRARY(PTHREAD_LIBRARY NAMES pthread)
#   ENDIF(NOT INTEL_GUIDE_LIBRARY MATCHES "NOTFOUND")
#   IF(NOT BLAS_guide_LIBRARY MATCHES "NOTFOUND")
#     SET(REQUIRE_PTHREAD TRUE)
#     FIND_LIBRARY(PTHREAD_LIBRARY NAMES pthread)
#   ENDIF(NOT BLAS_guide_LIBRARY MATCHES "NOTFOUND")
# 
#   IF(LAPACK_LIBRARY MATCHES "mkl")
#     MESSAGE(STATUS "Found intel/mkl library")
#     SET(LAPACK_LIBRARY_INIT 1)
#     SET(BLAS_LIBRARY_INIT 1)
#     SET(HAVE_MKL 1) # CACHE BOOL "HAVE_MKL is set to 1")
#     SET(MKL_INC_PATHS $ENV{MKL_HOME}/include ${MKL_PATHS}) 
#     FIND_PATH(MKL_INCLUDE_DIR mkl.h ${MKL_INC_PATHS})
#     INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
#   ENDIF(LAPACK_LIBRARY MATCHES "mkl")
# 
# ENDIF(NOT LAPACK_LIBRARY_INIT)


###################################################
# Looking for Veclib.
# (in case MKL was not found)
###################################################

IF(NOT LAPACK_LIBRARY_INIT)
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -framework vecLib")
    SET(LAPACK_LIBRARY_INIT 1)
    SET(MAC_VECLIB 1 CACHE BOOL "use Mac Framework")
    SET(LAPACK_LIBRARY "")
    MESSAGE(STATUS "Using Framework on Darwin.")

    ## check goto library: does not work so well
    #FIND_LIBRARY(BLAS_LIBRARY NAMES goto
    #    PATHS 
    #    $ENV{GOTOBLAS_HOME}
    #    /usr/lib
    #    /usr/local/lib
    #    /sw/lib
    #    )
    SET(BLAS_LIBRARY "")
    SET(BLAS_LIBRARY_INIT 1)
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
ENDIF(NOT LAPACK_LIBRARY_INIT)


###################################################
# Looking for ATLAS.
###################################################

IF(NOT LAPACK_LIBRARY_INIT)
  IF($ENV{ATLAS} MATCHES "atlas")
    IF($ENV{ATLAS} MATCHES "lapack") 
      SET(LAPACK_LIBRARY_INIT 1)
    ENDIF($ENV{ATLAS} MATCHES "lapack") 
    SET(BLAS_LIBRARY $ENV{ATLAS})
    SET(BLAS_LIBRARY_INIT 1)
  ENDIF($ENV{ATLAS} MATCHES "atlas")

  IF($ENV{LAPACK} MATCHES "lapack")
    SET(LAPACK_LIBRARY $ENV{LAPACK})
    SET(LAPACK_LIBRARY_INIT 1)
  ENDIF($ENV{LAPACK} MATCHES "lapack")

  IF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")
    SET(ELIB essl)
    IF(ENABLE_OMP)
      SET(ELIB esslsmp)
    ENDIF(ENABLE_OMP)
   
    IF(NOT LAPACK_LIBRARY_INIT)
      SET(LLIB lapack-SP4_${QMC_BITS} lapack)
      FIND_LIBRARY(LAPACK_LIBRARY  
        NAMES ${LLIB}
        PATHS /usr/apps/math/lapack/LAPACK
        lib
        )
      FIND_LIBRARY(BLAS_LIBRARY ${ELIB}
                   /usr/lib
                  )
    ENDIF(NOT LAPACK_LIBRARY_INIT)

    SET(LAPACK_LIBRARY_INIT 1)
    SET(BLAS_LIBRARY_INIT 1)
    MESSAGE(STATUS "Found lapack/blas on AIX system")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")
ENDIF(NOT LAPACK_LIBRARY_INIT)


###################################################
# Looking for CLAPACK on Windows.
###################################################

IF(NOT LAPACK_LIBRARY_INIT)
  if(WIN32 AND NOT UNIX)
    set(TRIAL_PATHS "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/lib")
    find_library(LAPACK_LIBRARY lapack ${TRIAL_PATHS})
    find_library(BLAS_LIBRARY blas ${TRIAL_PATHS})
    if(LAPACK_LIBRARY AND BLAS_LIBRARY)
      SET(LAPACK_LIBRARY_INIT 1)
      SET(BLAS_LIBRARY_INIT 1)
      set(LAPACK_DEFINITIONS "-DBOOST_NUMERIC_BINDINGS_USE_CLAPACK=1" CACHE STRING "Lapack definitions")
      MARK_AS_ADVANCED(LAPACK_DEFINITIONS)
      MESSAGE(STATUS "Found Lapack: ${LAPACK_LIBRARY} ${BLAS_LIBRARY}")
    endif(LAPACK_LIBRARY AND BLAS_LIBRARY)
  endif(WIN32 AND NOT UNIX)
ENDIF(NOT LAPACK_LIBRARY_INIT)


###################################################
# Falling back to CMake FindBLAS / LAPACK.
# (Fortran compiler is needed)
###################################################

IF(NOT BLAS_LIBRARY_INIT AND NOT LAPACK_LIBRARY_INIT)
  message(STATUS "Falling back to CMake provied LAPACK/BLAS detection.")
  find_package(BLAS)
  if(BLAS_FOUND)
    SET(BLAS_LIBRARY_INIT 1)
    SET(BLAS_LIBRARY ${BLAS_LIBRARIES})
    find_package(LAPACK)
    if(LAPACK_FOUND)
      SET(LAPACK_LIBRARY_INIT 1)
      SET(LAPACK_LIBRARY ${LAPACK_LIBRARIES})
    endif(LAPACK_FOUND)
  endif(BLAS_FOUND)
  if(NOT BLAS_LIBRARY_INIT AND NOT LAPACK_LIBRARY_INIT)
    message(STATUS "Enabling Fortran for LAPACK/BLAS detection.")
    enable_language(Fortran)
    find_package(BLAS)
    if(BLAS_FOUND)
      SET(BLAS_LIBRARY_INIT 1)
      SET(BLAS_LIBRARIES ${BLAS_LIBRARIES} -lgfortran)
      SET(BLAS_LIBRARY ${BLAS_LIBRARIES})
      find_package(LAPACK)
      if(LAPACK_FOUND)
        SET(LAPACK_LIBRARY_INIT 1)
        SET(LAPACK_LIBRARY ${LAPACK_LIBRARIES})
      endif(LAPACK_FOUND)
    endif(BLAS_FOUND)
  endif(NOT BLAS_LIBRARY_INIT AND NOT LAPACK_LIBRARY_INIT)
ENDIF(NOT BLAS_LIBRARY_INIT AND NOT LAPACK_LIBRARY_INIT)


###################################################
# Looking for SCALAPACK.
###################################################

IF(USE_SCALAPACK)
  SET(PNPATHS 
    ${MKL_PATHS}
    ${BLACS_HOME}/lib
    ${SCALAPACK_HOME}/lib
    /usr/lib
    /opt/lib
    /usr/local/lib
    /sw/lib
    )

  IF(INTEL_MKL)
    FIND_LIBRARY(BLACSLIB mkl_blacs_${PLAT}_lp${QMC_BITS} PATHS  ${PNPATHS})
    FIND_LIBRARY(SCALAPACKLIB mkl_scalapack PATHS  ${PNPATHS})
  ENDIF(INTEL_MKL)

  IF(NOT SCALAPACKLIB)
    FIND_LIBRARY(BLACSLIB blacs_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
    FIND_LIBRARY(BLACSCINIT blacsCinit_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
    FIND_LIBRARY(SCALAPACKLIB scalapack PATHS  ${PNPATHS})
  ENDIF(NOT SCALAPACKLIB)

  IF(BLACSLIB AND SCALAPACKLIB)
    SET(FOUND_SCALAPACK 1 CACHE BOOL "Found scalapack library")
  ELSE(BLACSLIB AND SCALAPACKLIB)
    SET(FOUND_SCALAPACK 0 CACHE BOOL "Mising scalapack library")
  ENDIF(BLACSLIB AND SCALAPACKLIB)

  MARK_AS_ADVANCED(
    BLACSCINIT
    BLACSLIB
    SCALAPACKLIB
    FOUND_SCALAPACK
    )
ENDIF(USE_SCALAPACK)



###################################################
# Finalize (setting some variables).
###################################################

message(STATUS "LAPACK DEBUG::LAPACK_LIBRARY = ${LAPACK_LIBRARY}")

if(BLAS_LIBRARY)
  SET(BLAS_LIBRARIES ${BLAS_LIBRARY})
endif(BLAS_LIBRARY)
if(LAPACK_LIBRARY)
  SET(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
endif(LAPACK_LIBRARY)

MARK_AS_ADVANCED(
  LAPACK_LIBRARY 
  BLAS_LIBRARY 
)

if(BLAS_LIBRARY_INIT)
  set(BLAS_FOUND TRUE)
endif(BLAS_LIBRARY_INIT)
if(LAPACK_LIBRARY_INIT)
  set(LAPACK_FOUND TRUE)
endif(LAPACK_LIBRARY_INIT)

IF(ALPS_BUILD_ON_CRAY)
  MARK_AS_ADVANCED(
    BLAS_FOUND
    LAPACK_FOUND
  )
  SET(BLAS_FOUND TRUE)
  SET(LAPACK_FOUND TRUE)
ENDIF(ALPS_BUILD_ON_CRAY)
