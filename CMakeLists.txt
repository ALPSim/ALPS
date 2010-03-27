#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include (CMakeDetermineSystem)
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ALPS_BUILD_EXAMPLES "Build ALPS examples" ON)
option(ALPS_BUILD_TESTS "Build ALPS tests" ON)
option(ALPS_BUILD_PYTHON "Build ALPS python extentions" ON)
option(ALPS_BUILD_APPLICATIONS "Build ALPS applications" ON)
option(ALPS_INSTALL_HEADERS "Install headers for ALPS and all dependent libraries" ON)
option(ALPS_INSTALL_DEPRECATED "Install deprecated headers for ALPS" OFF)
option(ALPS_BUILD_DEVELOPER_TOOLS "Build toold used by developers to maintain ALPS" OFF)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_BUILD_ON_CRAY "Build ALPS for cray" OFF)
option(ALPS_BUILD_DOC "Build ALPS documentation" OFF)


# a workaround for Mac OS X's LLVM compiler
option(ALPS_LLVM_WORKAROUND "Exclude some parts of ALPS from the build with LLVM" OFF)
mark_as_advanced(ALPS_LLVM_WORKAROUND)

#maintenance building options
option(ALPS_FOR_VISTRAILS "Build for use within Vistrails" OFF)
option(ALPS_PACKAGE_LIBRARIES "Install all dependent libraries into package" OFF)
option(ALPS_DOWNLOAD_PACKAGES "Prepare an installer that downloads instead of containing the installation packages" OFF)
option(ALPS_PATCH_VISTRAILS "Patch Vistrails from Vistrails SVN" OFF)

mark_as_advanced(ALPS_PATCH_VISTRAILS)

if (ALPS_FOR_VISTRAILS)
  set(ALPS_INSTALL_HEADERS OFF)
  set(ALPS_INSTALL_DEPRECATED OFF)
  set(ALPS_BUILD_DOC OFF)
else (ALPS_FOR_VISTRAILS)
  set(ALPS_PATCH_VISTRAILS OFF)
endif (ALPS_FOR_VISTRAILS)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################

if(WIN32 AND NOT UNIX)
  if (ALPS_FOR_VISTRAILS)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files (x86)/Vistrails" CACHE PATH "ALPS install prefix")
    set(CMAKE_INSTALL_ALTERNATE_PREFIX "$ENV{HOMEDRIVE}/Program Files/Vistrails" CACHE PATH "ALPS install prefix")
    set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the Vistrails application")
    set(VISTRAILS_PACKAGE_DIR packages)
    set(VISTRAILS_PYTHON_EXTENSION_DIR Python25/Lib/site-packages )
    set(VISTRAILS_LIB_DIR lib)
#    set(VISTRAILS_PYTHON_EXTENSION_DIR Python25/DLLs)
    set(VISTRAILS_DYLIB_DIR "")
  else (ALPS_FOR_VISTRAILS)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files (x86)/ALPS" CACHE PATH "ALPS install prefix" FORCE)
      set(CMAKE_INSTALL_ALTERNATE_PREFIX "$ENV{HOMEDRIVE}/Program Files/ALPS" CACHE PATH "ALPS install prefix")
    endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  endif (ALPS_FOR_VISTRAILS)
else(WIN32 AND NOT UNIX)
  if (ALPS_FOR_VISTRAILS)
    if (APPLE)
      if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "/Applications/Vistrails/" CACHE PATH "ALPS install prefix" FORCE)
      endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(VISTRAILS_APP_NAME "Vistrails.app" CACHE STRING "Name of the Vistrails application")
      set(VISTRAILS_PACKAGE_DIR Contents/Resources/lib/python2.5/packages)
      set(VISTRAILS_PYTHON_EXTENSION_DIR Contents/Resources/lib/python2.5)
      set(VISTRAILS_LIB_DIR Contents/Resources/lib)
      set(VISTRAILS_DYLIB_DIR Contents/Frameworks)
      set(CMAKE_OSX_DEPLOYMENT_TARGET "10.5")
    else (APPLE)
      set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the Vistrails application")
      set(VISTRAILS_PACKAGE_DIR vistrails/packages)
      set(VISTRAILS_PYTHON_EXTENSION_DIR /usr/lib/python2.5/site-packages)
      set(VISTRAILS_LIB_DIR vistrails/lib)
      set(VISTRAILS_DYLIB_DIR lib)
    endif (APPLE)
  else (ALPS_FOR_VISTRAILS)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(CMAKE_INSTALL_PREFIX "/opt/alps" CACHE PATH "ALPS install prefix" FORCE)
    endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  endif (ALPS_FOR_VISTRAILS)
endif(WIN32 AND NOT UNIX)

######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/CMake)

######################################################################
# Version information
######################################################################
set(ALPS_VERSION_MAJOR 2)
set(ALPS_VERSION_MINOR 0)
set(ALPS_VERSION_BUILD 0d)
set(ALPS_YEAR 2010)
set(ALPS_RELEASE_PATCH "")

IF (EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
  # Probably a SVN workspace, determine revision level 
  FIND_PACKAGE(Subversion) 
  IF(Subversion_FOUND) 
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} ALPS) 
    MESSAGE("Current revision is ${ALPS_WC_REVISION}") 
  ELSE(Subversion_FOUND) 
    IF (NOT ALPS_RELEASE_PATCH)
      SET(ALPS_WC_REVISION "4050")
    ELSE (NOT ALPS_RELEASE_PATCH)
      SET(ALPS_WC_REVISION ${ALPS_RELEASE_PATCH})
    ENDIF (NOT ALPS_RELEASE_PATCH)
  ENDIF(Subversion_FOUND) 
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

set(VERSION_BUILD ${ALPS_WC_REVISION})

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_BUILD}-r${VERSION_BUILD}")
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")

set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/lib/xml")
if (WIN32 AND NOT UNIX)
  set(ALPS_XML_ALTERNATE_DIR "${CMAKE_INSTALL_ALTERNATE_PREFIX}/lib/xml")
endif (WIN32 AND NOT UNIX)
set(ALPS_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(ALPS_SRCDIR "${CMAKE_SOURCE_DIR}")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(bindir "${CMAKE_INSTALL_PREFIX}/bin")

######################################################################
# C & C++ Headers
######################################################################
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(rpc/rpc.h ALPS_HAVE_RPC_XDR_H)
CHECK_INCLUDE_FILE(stdarg.h ALPS_HAVE_STDARG_H)
CHECK_INCLUDE_FILE(sys/stat.h ALPS_HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE(sys/systeminfo.h ALPS_HAVE_SYS_SYSTEMINFO_H)
CHECK_INCLUDE_FILE(sys/types.h ALPS_HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(unistd.h ALPS_HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h ALPS_HAVE_WINDOWS_H)

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(valarray ALPS_HAVE_VALARRAY)

######################################################################
# Packages
######################################################################
#set(PARCOLLET_INTERPRETER "python" CACHE PATH "Python interpreter")


#find_package(Parcollet)
#include_directories(${PARCOLLET_INCLUDE_DIR})
#link_libraries(${PARCOLLET_LIBRARIES})


# Boost Source File
find_package(BoostSrc REQUIRED)

# MPI
IF(NOT MPI_FOUND)
  find_package(MPI)
ENDIF(NOT MPI_FOUND)
set(MPI_DEFINITIONS ${MPI_COMPILE_FLAGS})
if(MPI_INCLUDE_DIR)
  set(MPI_DEFINITIONS "${MPI_COMPILE_FLAGS} -I${MPI_INCLUDE_DIR}")
endif(MPI_INCLUDE_DIR)
  
set(MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH})
IF(MPI_FOUND AND NOT (ALPS_FOR_VISTRAILS AND WIN32 AND NOT LINUX))
  LINK_LIBRARIES(${MPI_LIBRARIES})
  include_directories(${MPI_INCLUDE_DIR})
  SET(ALPS_HAVE_MPI 1)
  SET(ALPS_HAVE_BOOST_MPI 1)
  set(BUILD_BOOST_MPI TRUE)
ELSE(MPI_FOUND AND NOT (ALPS_FOR_VISTRAILS AND WIN32 AND NOT LINUX))
  set(BUILD_BOOST_MPI FALSE)
ENDIF(MPI_FOUND AND NOT (ALPS_FOR_VISTRAILS AND WIN32 AND NOT LINUX))
IF(MPI_FOUND)
  foreach(arg ${MPI_LIBRARIES})
    set(MPI_LIBS "${MPI_LIBS} ${arg}")
endforeach(arg ${MPI_LIBRARIES})
ENDIF(MPI_FOUND)
#add_definitions( -DBOOST_MPI_HOMOGENEOUS )

# OpenMP
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(ALPS_ENABLE_OPENMP)

# Lapack
if (NOT ALPS_FOR_VISTRAILS)
  if (NOT LAPACK_FOUND)
    find_package(Lapack)
  endif (NOT LAPACK_FOUND)
  if (LAPACK_FOUND)
    set(ALPS_HAVE_LAPACK 1)
    set(ALPS_HAVE_BLAS 1)
  endif(LAPACK_FOUND)
  if(MSVC)
    file(COPY ${LAPACK_LIBRARY} ${BLAS_LIBRARY} DESTINATION ${PROJECT_BINARY_DIR}/lib)
  endif(MSVC)
  IF(NOT WIN32)
    IF (REQUIRE_PTHREAD)
      SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
    ENDIF (REQUIRE_PTHREAD)
  ENDIF(NOT WIN32)
endif (NOT ALPS_FOR_VISTRAILS)
if(MAC_VECLIB)
  set(LAPACK_LDFLAGS "-faltivec -framework vecLib")
endif(MAC_VECLIB) 
set(LAPACK_LINKER_FLAGS ${LAPACK_LDFLAGS})

# SQLite
if (NOT ALPS_FOR_VISTRAILS)
  find_package(SQLite)
  if(MSVC)
    file(COPY ${SQLite_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
    file(COPY ${SQLite_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
endif (NOT ALPS_FOR_VISTRAILS)

# LPSolve
if (NOT ALPS_FOR_VISTRAILS)
  find_package(LPSolve)
  if(MSVC)
    file(COPY ${LPSolve_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
    file(COPY ${LPSolve_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
endif (NOT ALPS_FOR_VISTRAILS)

# Zlib
find_package(ZLIB)
if(ZLIB_FOUND AND NOT ALPS_FOR_VISTRAILS)
  include_directories(${ZLIB_INCLUDE_DIR})
  link_libraries(${ZLIB_LIBRARIES})
  if(MSVC)
    file(COPY ${ZLIB_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
    file(COPY ${ZLIB_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
endif(ZLIB_FOUND AND NOT ALPS_FOR_VISTRAILS)
  
# SZip
find_package(SZIP)
if(SZIP_FOUND AND NOT ALPS_FOR_VISTRAILS)
  include_directories(${SZIP_INCLUDE_DIR})
  link_libraries(${SZIP_LIBRARIES})
  if(MSVC)
    file(COPY ${SZIP_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
    file(COPY ${SZIP_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
endif(SZIP_FOUND AND NOT ALPS_FOR_VISTRAILS)
  
# HDF5
find_package(HDF5)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)  
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARIES})
  SET(ALPS_HAVE_HDF5 1)
  IF (NOT ALPS_FOR_VISTRAILS)
    if(MSVC)
      file(COPY ${HDF5_hdf5_LIBRARY_RELEASE} ${HDF5_hdf5_hl_LIBRARY_RELEASE}
        DESTINATION ${PROJECT_BINARY_DIR}/lib)
      file(COPY ${HDF5_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
    endif(MSVC)
  ENDIF(NOT ALPS_FOR_VISTRAILS)
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "Require hdf5 1.8.2 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)
IF(HDF5_IS_PARALLEL)
  SET(ALPS_HAVE_HDF5_PARALLEL 1)
  MESSAGE(WARNING "parallel(MPI) hdf5 is detected. We will compile but ALPS does not use parallel HDF5. The standard version is preferred.")
  IF(NOT MPI_FOUND)
    MESSAGE(FATAL_ERROR "parallel(MPI) hdf5 needs MPI. Enable MPI or install serial HDF5 libraries.")
  ENDIF(NOT MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(HDF5_IS_PARALLEL)

# python
if(WIN32 AND NOT UNIX AND ALPS_BUILD_PYTHON)
  find_package(PyInstaller)
  include(UsePyInstaller)
endif(WIN32 AND NOT UNIX AND ALPS_BUILD_PYTHON)
find_package(PythonInterp)
set(PYTHON_SCRIPTDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
find_package(PythonLibs)


IF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
  find_package(Numpy)
    if (ALPS_FOR_VISTRAILS)
      if(UNIX AND NOT APPLE AND NOT WIN32)
        SET(PYTHON_INCLUDE_DIRS /usr/include/python2.5)
        SET(PYTHON_LIBRARY /usr/lib/libpython2.5.so)
      ELSE(UNIX AND NOT APPLE AND NOT WIN32)
        if(WIN32 AND NOT UNIX)
          SET(PYTHON_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Python25/include)
          SET(PYTHON_LIBRARY "${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Python25/libs/python25.lib")
        ELSE(WIN32 AND NOT UNIX)
          SET(PYTHON_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Contents/Frameworks/Python.framework/Headers)
          SET(PYTHON_LIBRARY "-F${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Contents/Frameworks -framework Python")
        endif(WIN32 AND NOT UNIX)
      endif(UNIX AND NOT APPLE AND NOT WIN32)
      set(PYTHONLIBS_FOUND TRUE)
    endif(ALPS_FOR_VISTRAILS)
  include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
  SET(ALPS_HAVE_PYTHON 1)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS}) 
  link_libraries(${PYTHON_LIBRARY})
ELSE (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
  set(BUILD_BOOST_PYTHON FALSE)
ENDIF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)

######################################################################
# configuration files
######################################################################
configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/src/alps/config.h)
configure_file(src/alps/version.h.in ${PROJECT_BINARY_DIR}/src/alps/version.h)

configure_file(lib/xml/lattices.xml.in ${PROJECT_BINARY_DIR}/lib/xml/lattices.xml)
configure_file(lib/xml/models.xml.in ${PROJECT_BINARY_DIR}/lib/xml/models.xml)

######################################################################
# include directories
######################################################################
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(src)
include_directories(${Boost_ROOT_DIR})

######################################################################
# RPATH setting
######################################################################
if(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

######################################################################
# install dependent libraries under Windows, from applications
######################################################################

if(WIN32 AND NOT UNIX AND NOT ALPS_FOR_VISTRAILS)
  # Lapack
  if(LAPACK_LIBRARIES)
    foreach(file ${LAPACK_LIBRARIES})
      get_filename_component(LAPACK_PATH ${file} PATH)
    endforeach(file ${LAPACK_LIBRARIES})
    add_custom_target(liblapack ALL
      COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LAPACK_PATH}/lapack.lib ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LAPACK_PATH}/blas.lib ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LAPACK_PATH}/f2c.lib ${PROJECT_BINARY_DIR}/lib
    )
    install(FILES ${PROJECT_BINARY_DIR}/lib/lapack.lib ${PROJECT_BINARY_DIR}/lib/blas.lib
            ${PROJECT_BINARY_DIR}/lib/f2c.lib
            DESTINATION lib COMPONENT dependentlibraries
    )
  endif(LAPACK_LIBRARIES)
  # SQLite
  if(SQLite_LIBRARIES)
    get_filename_component(SQLite_PATH ${SQLite_LIBRARIES} PATH)
    get_filename_component(SQLite_BASE ${SQLite_PATH} PATH)
    get_filename_component(SQLite_NAME ${SQLite_LIBRARIES} NAME_WE)
    add_custom_target(libsqlite ALL
      COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${SQLite_LIBRARIES} ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${SQLite_BASE}/bin/${SQLite_NAME}.dll ${PROJECT_BINARY_DIR}/bin
    )
    install(FILES ${PROJECT_BINARY_DIR}/lib/${SQLite_NAME}.lib DESTINATION lib COMPONENT dependentlibraries)
    install(FILES ${PROJECT_BINARY_DIR}/bin/${SQLite_NAME}.dll DESTINATION bin COMPONENT dependentlibraries)
  endif(SQLite_LIBRARIES)
  # LPSolve
  if(LPSolve_LIBRARIES)
    get_filename_component(LPSolve_PATH ${LPSolve_LIBRARIES} PATH)
    get_filename_component(LPSolve_BASE ${LPSolve_PATH} PATH)
    get_filename_component(LPSolve_NAME ${LPSolve_LIBRARIES} NAME_WE)
    add_custom_target(liblpsolve ALL
      COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LPSolve_LIBRARIES} ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LPSolve_BASE}/bin/${LPSolve_NAME}.dll ${PROJECT_BINARY_DIR}/bin
    )
    install(FILES ${PROJECT_BINARY_DIR}/lib/${LPSolve_NAME}.lib DESTINATION lib COMPONENT dependentlibraries)
    install(FILES ${PROJECT_BINARY_DIR}/bin/${LPSolve_NAME}.dll DESTINATION bin COMPONENT dependentlibraries)
  endif(LPSolve_LIBRARIES)
  # alps & boost
  add_custom_target(libalps ALL
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib/xml
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/archive.exe ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/parameter2xml.exe ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/ising.exe ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E copy_directory ${CMAKE_INSTALL_PREFIX}/lib/xml ${PROJECT_BINARY_DIR}/lib/xml
  )
endif(WIN32 AND NOT UNIX AND NOT ALPS_FOR_VISTRAILS)

######################################################################
# setting up ALPSConfig.cmake and include.mk
######################################################################

set(ALPS_USE_FILE ${CMAKE_INSTALL_PREFIX}/share/alps/UseALPS.cmake)

set(ALPS_INCLUDE_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/include)
# if(NOT INSTALL_BOOST_HEADERS)
  list(APPEND ALPS_INCLUDE_DIRS_CONFIG ${Boost_ROOT_DIR})
# endif(NOT INSTALL_BOOST_HEADERS)
set(ALPS_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)
set(ALPS_RUNTIME_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/bin)

if(MPI_DEFINITIONS)
  set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${MPI_DEFINITIONS}")
endif(MPI_DEFINITIONS)
if(LAPACK_DEFINITIONS)
  set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${LAPACK_DEFINITIONS}")
endif(LAPACK_DEFINITIONS)

if(MPI_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
endif(MPI_INCLUDE_DIR)
if(SQLite_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${SQLite_INCLUDE_DIR})
endif(SQLite_INCLUDE_DIR)
if(LPSolve_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${LPSolve_INCLUDE_DIR})
endif(LPSolve_INCLUDE_DIR)
if(HDF5_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
endif(HDF5_INCLUDE_DIR)

if(MPI_LINKER_FLAGS)
  set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
endif(MPI_LINKER_FLAGS)
if(LAPACK_LINKER_FLAGS)
  set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
endif(LAPACK_LINKER_FLAGS)

if(MPI_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${MPI_LIBRARIES})
endif(MPI_LIBRARIES)
if(LAPACK_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${LAPACK_LIBRARIES})
endif(LAPACK_LIBRARIES)
if(SQLite_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${SQLite_LIBRARIES})
endif(SQLite_LIBRARIES)
if(LPSolve_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${LPSolve_LIBRARIES})
endif(LPSolve_LIBRARIES)
if(HDF5_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
endif(HDF5_LIBRARIES)

set (HDF5_LIBS "")
foreach (l ${HDF5_LIBRARIES})
  set (HDF5_LIBS "${HDF5_LIBS} ${l}")
  if (${l} STREQUAL "debug")
    set (HDF5_LIBS "")
  endif()
  if (${l} STREQUAL "optimized")
    set (HDF5_LIBS "")
  endif()
endforeach(l)

set (LAPACK_LIBS "")
foreach (l ${LAPACK_LIBRARY})
  set (LAPACK_LIBS "${LAPACK_LIBS} ${l}")
endforeach(l)

configure_file(config/CMake/ALPSConfig.cmake.in ${PROJECT_BINARY_DIR}/config/CMake/ALPSConfig.cmake)
configure_file(config/CMake/include.mk.in ${PROJECT_BINARY_DIR}/config/CMake/include.mk)

set (HDF5_LIBRARIES ${HDF5_LIBRARIES_SAVE})


# installation
######################################################################
include(InstallRequiredSystemLibraries)

if (ALPS_INSTALL_HEADERS)
  set(ALPS_HEADER_DIR "include")
  install(DIRECTORY src/alps src/boost src/ietl src/mocasito COMPONENT headers DESTINATION ${ALPS_HEADER_DIR}
          FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h ${PROJECT_BINARY_DIR}/src/alps/version.h DESTINATION ${ALPS_HEADER_DIR}/alps COMPONENT headers)
  if (ALPS_INSTALL_DEPRECATED)
    install(DIRECTORY legacy/alps COMPONENT deprecated DESTINATION ${ALPS_HEADER_DIR}
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
            PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  endif (ALPS_INSTALL_DEPRECATED)

  if(ALPS_PACKAGE_LIBRARIES)
     install(DIRECTORY ${Boost_ROOT_DIR}/boost DESTINATION ${ALPS_HEADER_DIR} COMPONENT dependentheaders
             FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
             PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
     install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/include DESTINATION . COMPONENT dependentheaders)
  endif(ALPS_PACKAGE_LIBRARIES)

endif(ALPS_INSTALL_HEADERS)



if (ALPS_FOR_VISTRAILS)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR})
else (ALPS_FOR_VISTRAILS)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION lib COMPONENT xml)
endif (ALPS_FOR_VISTRAILS)

if(ALPS_PACKAGE_LIBRARIES)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/lib DESTINATION . COMPONENT dependentlibraries)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/bin DESTINATION . COMPONENT dependentlibraries)
endif(ALPS_PACKAGE_LIBRARIES)

if (ALPS_FOR_VISTRAILS)
  if (NOT WIN32 AND NOT APPLE)
    install(DIRECTORY lib/pyalps DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}
          FILES_MATCHING PATTERN "*.py"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)
  else (NOT WIN32 AND NOT APPLE)
    install(DIRECTORY lib/pyalps DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}
            FILES_MATCHING PATTERN "*.py"
            PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)
  endif (NOT WIN32 AND NOT APPLE)
  install(DIRECTORY vistrails/alps DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PACKAGE_DIR}
        FILES_MATCHING PATTERN "*.py"  PATTERN "*.xml"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  install(DIRECTORY lib/xml DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR}
        FILES_MATCHING PATTERN "*.xsl"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
else (ALPS_FOR_VISTRAILS)
  install(DIRECTORY lib/xml DESTINATION lib COMPONENT xml
        FILES_MATCHING PATTERN "*.xsl"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION lib COMPONENT xml)

  install(FILES config/CMake/UseALPS.cmake 
          ${PROJECT_BINARY_DIR}/config/CMake/ALPSConfig.cmake
          ${PROJECT_BINARY_DIR}/config/CMake/include.mk
          DESTINATION share/alps COMPONENT build)

  install(DIRECTORY lib/pyalps DESTINATION lib  COMPONENT python
        FILES_MATCHING PATTERN "*.py"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
endif (ALPS_FOR_VISTRAILS)

######################################################################
# patch Vistrails
######################################################################

if(ALPS_PATCH_VISTRAILS)
  find_package(VistrailsSVN)
  MESSAGE (STATUS "patching Vistrails from ${VISTRAILS_SVN_DIR}")
  if (APPLE)
    foreach (DIR api core db gui packages tests)
      install(DIRECTORY ${VISTRAILS_SVN_DIR}/vistrails/${DIR} 
              DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}
              PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.pyc" EXCLUDE)
    endforeach(DIR)
    install(FILES ${VISTRAILS_SVN_DIR}/vistrails/vistrails.py 
            DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR/../..})
  else(APPLE)
    if (WIN32) 
        install(DIRECTORY ${VISTRAILS_SVN_DIR}/vistrails
                DESTINATION .
                PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.pyc" EXCLUDE)
    else (WIN32)
      message(ERROR "Cannot patch VisTrails on your systems yet. Please add the needed code.")
    endif (WIN32)
  endif(APPLE)
endif(ALPS_PATCH_VISTRAILS)


######################################################################
# documentation
######################################################################

if (ALPS_BUILD_DOC AND NOT ALPS_FOR_VISTRAILS AND NOT WIN32)
add_subdirectory(doc)
endif (ALPS_BUILD_DOC AND NOT ALPS_FOR_VISTRAILS AND NOT WIN32)

######################################################################
# libraries
######################################################################

add_subdirectory(src/alps)
add_subdirectory(src/boost)
add_subdirectory(lib/pyalps)

######################################################################
# programs
######################################################################

add_subdirectory(tool)
if (NOT ALPS_FOR_VISTRAILS)
  add_subdirectory(example)
endif (NOT ALPS_FOR_VISTRAILS)

######################################################################
# applications
######################################################################

if (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_APPLICATIONS)

######################################################################
# developer tools
######################################################################

if (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_DEVELOPER_TOOLS)
  add_subdirectory(script)
endif (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_DEVELOPER_TOOLS)
    

##########################################################################
# Binary packages                                                        #
##########################################################################

 set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "ALPS and Boost libraries")
 set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "ALPS headers")
 set(CPACK_COMPONENT_DEPRECATED_DISPLAY_NAME "Deprecated ALPS headers")
 set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "ALPS applications")
 set(CPACK_COMPONENT_XML_DISPLAY_NAME "XML files")
 set(CPACK_COMPONENT_PYTHON_DISPLAY_NAME "ALPS Python packages")
 set(CPACK_COMPONENT_TOOLS_DISPLAY_NAME "ALPS tools")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_DISPLAY_NAME "Dependent libraries")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_DISPLAY_NAME "Headers of dependent libraries")
 set(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "ALPS examples")
 set(CPACK_COMPONENT_BUILD_DISPLAY_NAME "Build system support")
 set(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
 
 set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The libraries needed by the ALPS applications")
 set(CPACK_COMPONENT_HEADERS_DESCRIPTION "The ALPS header files")
 set(CPACK_COMPONENT_DEPRECATED_DESCRIPTION "Deprecated ALPS header files")
 set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION "The ALPS application programs")
 set(CPACK_COMPONENT_XML_DESCRIPTION "The XML files and XSL stylefiles for ALPS input and evaluation")
 set(CPACK_COMPONENT_PYTHON_DESCRIPTION "The ALPS Python packages for data evaluation")
 set(CPACK_COMPONENT_TOOLS_DESCRIPTION "The ALPS tools")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_DESCRIPTION "External libraries that ALPS depends on. They do not need to be installed if they already exist on the system and can be found by the ALPS applications")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_DESCRIPTION "Headers of libraries that ALPS depends on")
 set(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "ALPS examples")
 set(CPACK_COMPONENT_BUILD_DESCRIPTION "Include files for CMake and make to simplify the building of ALPS applications")
 set(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation for the ALPS libraries")

 set(CPACK_COMPONENT_LIBRARIES_GROUP "Runtime")
 set(CPACK_COMPONENT_APPLICATIONS_GROUP "Runtime")
 set(CPACK_COMPONENT_XML_GROUP "Runtime")
 set(CPACK_COMPONENT_PYTHON_GROUP "Runtime")
 set(CPACK_COMPONENT_TOOLS_GROUP "Runtime")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_GROUP "Runtime")
 set(CPACK_COMPONENT_HEADERS_GROUP "Development")
 set(CPACK_COMPONENT_DEPRECATED_GROUP "Development")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_GROUP "Development")
 set(CPACK_COMPONENT_BUILD_GROUP "Development")
 #set(CPACK_COMPONENT_EXAMPLES_GROUP "Examplecodes")

 set(CPACK_COMPONENT_HEADERS_DISABLED ON)
 set(CPACK_COMPONENT_DEPRECATED_DISABLED ON)
 set(CPACK_COMPONENT_DEPENDENTHEADERS_DISABLED ON)
 set(CPACK_COMPONENT_BUILD_DISABLED ON)

 set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "Components needed only to build applications based on the ALPS libraries")
 set(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "Components needed to run the ALPS applications")
# set(CPACK_COMPONENT_GROUP_EXAMPLECODES_DESCRIPTION "Installs some example programs")
 set(CPACK_COMPONENT_APPLICATIONS_DEPENDS libraries)
 set(CPACK_COMPONENT_PYTHON_DEPENDS libraries)
 set(CPACK_COMPONENT_TOOLS_DEPENDS libraries)
 set(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries)
 
 
if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_NAME "alps-vistrails-patch")
  set(CPACK_NSIS_MODIFY_PATH OFF)
else(ALPS_FOR_VISTRAILS)
  set(CPACK_NSIS_MODIFY_PATH ON)
  if(ALPS_DOWNLOAD_PACKAGES)
    set(CPACK_PACKAGE_NAME "alps-download")
  else(ALPS_DOWNLOAD_PACKAGES)
    set(CPACK_PACKAGE_NAME "alps")
  endif(ALPS_DOWNLOAD_PACKAGES)
endif(ALPS_FOR_VISTRAILS)
set(CPACK_PACKAGE_VENDOR "ALPS Collaboration")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ALPS ${ALPS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")

if (ALPS_PACKAGE_LIBRARIES)
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README-package.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE-package.txt")
else (ALPS_PACKAGE_LIBRARIES)
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/applications/LICENSE.txt")
endif (ALPS_PACKAGE_LIBRARIES)

set(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/Welcome.txt")
set(CPACK_PACKAGE_VERSION "${ALPS_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${ALPS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${ALPS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_SUBMINOR "${ALPS_VERSION_SUBMINOR}")

if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "VisTrails")
else (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "ALPS")
endif (ALPS_FOR_VISTRAILS)
set(CPACK_PACKAGING_INSTALL_PREFIX "/")

set(CPACK_PACKAGE_DEFAULT_LOCATION ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_VERSION_PATCH "${ALPS_VERSION_BUILD}-${VERSION_BUILD}")

# For Debian package
set(CPACK_DEBIAN_PACKAGE_NAME "alps")
set(CPACK_DEBIAN_PACKAGE_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}~d1")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
# set(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS)
# set(CPACK_DEBIAN_PACKAGE_DEPENDS)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>")
# set(CPACK_DEBIAN_PACKAGE_DESCRIPTION)
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
# set(CPACK_DEBIAN_PACKAGE_RECOMMENDS)
# set(CPACK_DEBIAN_PACKAGE_SUGGESTS)

# For RPM package
# set(CPACK_RPM_PACKAGE_SUMMARY "${CPACK_PACKAGE_DESCRIPTION}")
# set(CPACK_RPM_PACKAGE_NAME "${CPACK_PACKAGE_NAME}")
# set(CPACK_RPM_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "amd64")
# set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_LICENCE "ALPS Licence")
set(CPACK_RPM_PACKAGE_GROUP "${CPACK_PACKAGE_VENDOR}")
# set(CPACK_RPM_PACKAGE_VENDOR "${CPACK_PACKAGE_VENDOR}")
# set(CPACK_RPM_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION_FILE}")
# set(CPACK_RPM_PACKAGE_REQUIRES "")
# set(CPACK_RPM_SPEC_INSTALL_POST "")
# set(CPACK_RPM_SPEC_MORE_DEFINE "")
# set(CPACK_RPM_USER_SPECFILE "")
# set(CPACK_RPM_GENERATE_USER_BINARY_SPECFILE_TEMPLATE "")
# set(CPACK_RPM_PRE_INSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_PRE_UNINSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_POST_UNINSTALL_SCRIPT_FILE "")

if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  # NOTE: No ALPS icon yet
#  set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/tools/build/CMake\\\\InstallIcon.bmp")
#  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}")
  set(CPACK_NSIS_CONTACT "wistaria@comp-phys.org")
  set(CPACK_NSIS_HELP_LINK "http://alps.comp-phys.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http://alps.comp-phys.org")
  
  # Encode the compiler name in the package 
#  if(NOT ALPS_FOR_VISTRAILS)
#    if (MSVC90)
#      set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc9")
#      set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME}")
#    elseif (MSVC90)
#      set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-unknown")  
#      set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for unknown compiler")    
#    endif (MSVC90)
#  endif(NOT ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_NSIS_DISPLAY_NAME}")
endif(WIN32 AND NOT UNIX)

#set(CPACK_STRIP_FILES ON)
include(CPack)

if (ALPS_DOWNLOAD_PACKAGES)
  if(COMMAND cpack_configure_downloads)
    cpack_configure_downloads(
    "http://alps.comp-phys.org/static/software/releases/packages/"
      ALL ADD_REMOVE)
  endif()
endif(ALPS_DOWNLOAD_PACKAGES)

##########################################################################


######################################################################
# tests
######################################################################

if (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_TESTS)
  ENABLE_TESTING()
  add_subdirectory(test)
endif (NOT ALPS_FOR_VISTRAILS AND ALPS_BUILD_TESTS)
