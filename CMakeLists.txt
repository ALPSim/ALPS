cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# option(INSTALL_BOOST_HEADERS "Install all the Boost header files" OFF)

if(MSVC)
  SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /openmp" CACHE STRING
    "Flags used by the compiler." FORCE)
endif(MSVC)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(WIN32 AND NOT UNIX)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/ALPS" CACHE PATH "ALPS install prefix" FORCE)
  else(WIN32 AND NOT UNIX)
    set(CMAKE_INSTALL_PREFIX "/opt/alps" CACHE PATH "ALPS install prefix" FORCE)
  endif(WIN32 AND NOT UNIX)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/CMake)

######################################################################
# Version information
######################################################################
set(ALPS_VERSION_MAJOR 1)
set(ALPS_VERSION_MINOR 3)
set(ALPS_VERSION_BUILD 5)
set(ALPS_YEAR 2009)

set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_BUILD}")
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")

set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/lib/xml")

######################################################################
# C & C++ Headers
######################################################################
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(rpc/rpc.h ALPS_HAVE_RPC_XDR_H)
CHECK_INCLUDE_FILE(stdarg.h ALPS_HAVE_STDARG_H)
CHECK_INCLUDE_FILE(sys/systeminfo.h ALPS_HAVE_SYS_SYSTEMINFO_H)
CHECK_INCLUDE_FILE(sys/types.h ALPS_HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(unistd.h ALPS_HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h ALPS_HAVE_WINDOWS_H)

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(valarray ALPS_HAVE_VALARRAY)

######################################################################
# Packages
######################################################################

# Boost Source File
find_package(BoostSrc REQUIRED)

# MPI
find_package(MPI)
IF(MPI_FOUND)
  LINK_LIBRARIES(${MPI_LIBRARIES})
ENDIF(MPI_FOUND)

# Lapack
find_package(Lapack)

# SQLite
find_package(SQLite)

# LPSolve
find_package(LPSolve)

# Zlib
find_package(ZLIB)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIR})
  link_libraries(${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)
  
# SZip
find_package(SZIP)
if(SZIP_FOUND)
  include_directories(${SZIP_INCLUDE_DIR})
  link_libraries(${SZIP_LIBRARIES})
endif(SZIP_FOUND)
  
# HDF5
find_package(HDF5)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARIES})
  SET(ALPS_HAVE_HDF5 1)
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "Require hdf5 1.8.2 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)
IF(HDF5_IS_PARALLEL)
  SET(ALPS_HAVE_HDF5_PARALLEL 1)
  MESSAGE(WARNING "parallel(MPI) hdf5 is detected. We try to compile but ALPS do not use these feature.")
  IF(NOT MPI_FOUND)
    MESSAGE(FATAL_ERROR "parallel(MPI) hdf5 needs MPI. Enable MPI or install seriall hdf5.")
  ENDIF(NOT MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(HDF5_IS_PARALLEL)

# python
if(WIN32 AND NOT UNIX)
  find_package(PyInstaller)
  include(UsePyInstaller)
else(WIN32 AND NOT UNIX)
  find_package(PythonInterp)
  set(PYTHON_SCRIPTDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
endif(WIN32 AND NOT UNIX)

######################################################################
# configuration files
######################################################################
configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/src/alps/config.h)

configure_file(lib/xml/lattices.xml.in ${PROJECT_BINARY_DIR}/lib/xml/lattices.xml)
configure_file(lib/xml/models.xml.in ${PROJECT_BINARY_DIR}/lib/xml/models.xml)
configure_file(lib/xml/changestylesheet.xsl.in ${PROJECT_BINARY_DIR}/lib/xml/changestylesheet.xsl)
configure_file(lib/xml/plot2xsl.xsl.in ${PROJECT_BINARY_DIR}/lib/xml/plot2xsl.xsl)

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
# setting up ALPSConfig.cmake
######################################################################
set(ALPS_INCLUDE_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/include)
# if(NOT INSTALL_BOOST_HEADERS)
  list(APPEND ALPS_INCLUDE_DIRS_CONFIG ${Boost_ROOT_DIR})
# endif(NOT INSTALL_BOOST_HEADERS)

set(ALPS_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)
set(ALPS_RUNTIME_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/bin)

set(ALPS_USE_FILE ${CMAKE_INSTALL_PREFIX}/share/alps/UseALPS.cmake)

configure_file(config/CMake/ALPSConfig.cmake.in ${PROJECT_BINARY_DIR}/config/CMake/ALPSConfig.cmake)

######################################################################
# installation
######################################################################
include(InstallRequiredSystemLibraries)

set(ALPS_HEADER_DIR "include")
install(DIRECTORY src/alps src/boost DESTINATION ${ALPS_HEADER_DIR}
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h DESTINATION ${ALPS_HEADER_DIR}/alps)
# if(INSTALL_BOOST_HEADERS)
#   install(DIRECTORY ${Boost_ROOT_DIR}/boost DESTINATION ${ALPS_HEADER_DIR}
#           FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
#           PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
# endif(INSTALL_BOOST_HEADERS)

install(DIRECTORY lib/xml DESTINATION lib
        FILES_MATCHING PATTERN "*.xsl"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION lib)

install(FILES config/CMake/UseALPS.cmake ${PROJECT_BINARY_DIR}/config/CMake/ALPSConfig.cmake
        DESTINATION share/alps)

install(DIRECTORY lib/pyalps DESTINATION lib
        FILES_MATCHING PATTERN "*.py"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)

######################################################################
# libraries
######################################################################
add_subdirectory(src/alps)
add_subdirectory(src/boost)

######################################################################
# programs
######################################################################
add_subdirectory(tool)
# add_subdirectory(example/parapack/single)
# add_subdirectory(example/scheduler)

##########################################################################
# Binary packages                                                        #
##########################################################################
set(CPACK_PACKAGE_NAME "ALPS")
set(CPACK_PACKAGE_VENDOR "alps.comp-phys.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ALPS ${ALPS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/Welcome.txt")
set(CPACK_PACKAGE_VERSION "${ALPS_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${ALPS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${ALPS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_SUBMINOR "${ALPS_VERSION_SUBMINOR}")
set(CPACK_PACKAGING_INSTALL_PREFIX "/")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "/opt/alps")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/opt/alps")

if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  # NOTE: No ALPS icon yet
#  set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/tools/build/CMake\\\\InstallIcon.bmp")
#  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_SUBMINOR}")
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\alps.comp-phys.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\alps.comp-phys.org")
  set(CPACK_NSIS_CONTACT "comp-phys-alps-users@lists.comp-phys.org")
  set(CPACK_NSIS_MODIFY_PATH ON)
  
  # Encode the compiler name in the package 
  if (MSVC60)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc6")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual C++ 6")
  elseif (MSVC70)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc7")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2002")
  elseif (MSVC71)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc71")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2003")
  elseif (MSVC80)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc8")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2005")    
  elseif (MSVC90)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-vc9")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2008")
  elseif (BORLAND)
    set(CPACK_PACKAGE_FILE_NAME "ALPS-${ALPS_VERSION}-borland")  
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Borland C++ Builder")    
  endif (MSVC60)
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_NSIS_DISPLAY_NAME}")
endif(WIN32 AND NOT UNIX)
include(CPack)


if(COMMAND cpack_configure_downloads)
  cpack_configure_downloads(
  "http://alps.comp-phys.org/static/software/releases/"
    ALL ADD_REMOVE)
endif()
##########################################################################


######################################################################
# tests
######################################################################
ENABLE_TESTING()
add_subdirectory(test)
