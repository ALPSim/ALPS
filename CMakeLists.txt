include (CMakeDetermineSystem)
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ALPS_BUILD_EXAMPLES "Build ALPS examples" ON)
option(ALPS_BUILD_APPLICATIONS "Build ALPS applications" ON)
option(ALPS_INSTALL_HEADERS "Install headers for ALPS and all dependent libraries" ON)

#maintenance building options
option(ALPS_FOR_VISTRAILS "Build for use within Vistrails" OFF)
option(ALPS_PACKAGE_LIBRARIES "Install all dependent libraries into package" OFF)

if (ALPS_FOR_VISTRAILS)
  set(ALPS_INSTALL_HEADERS OFF)
endif (ALPS_FOR_VISTRAILS)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################

if(WIN32 AND NOT UNIX)
  if (ALPS_FOR_VISTRAILS)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/Vistrails" CACHE PATH "ALPS install prefix" FORCE)
    endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the Vistrails application")
    set(VISTRAILS_PACKAGE_DIR packages)
    set(VISTRAILS_PYTHON_EXTENSION_DIR Python25/Lib/site-packages )
    set(VISTRAILS_LIB_DIR lib)
#    set(VISTRAILS_PYTHON_EXTENSION_DIR Python25/DLLs)
    set(VISTRAILS_DYLIB_DIR "")
  else (ALPS_FOR_VISTRAILS)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/ALPS" CACHE PATH "ALPS install prefix" FORCE)
    endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  endif (ALPS_FOR_VISTRAILS)
else(WIN32 AND NOT UNIX)
  if (ALPS_FOR_VISTRAILS)
    message (status setting)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(CMAKE_INSTALL_PREFIX "/Applications/Vistrails/" CACHE PATH "ALPS install prefix" FORCE)
    endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(VISTRAILS_APP_NAME "Vistrails.app" CACHE STRING "Name of the Vistrails application")
    set(VISTRAILS_PACKAGE_DIR Contents/Resources/lib/python2.5/packages)
    set(VISTRAILS_PYTHON_EXTENSION_DIR Contents/Resources/lib/python2.5)
    set(VISTRAILS_LIB_DIR Contents/Resources/lib)
    set(VISTRAILS_DYLIB_DIR Contents/Frameworks)
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
set(ALPS_VERSION_MAJOR 1)
set(ALPS_VERSION_MINOR 4)
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
      SET(ALPS_WC_REVISION "svn")
    ELSE (NOT ALPS_RELEASE_PATCH)
      SET(ALPS_WC_REVISION ${ALPS_RELEASE_PATCH})
    ENDIF (NOT ALPS_RELEASE_PATCH)
  ENDIF(Subversion_FOUND) 
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

set(VERSION_BUILD ${ALPS_WC_REVISION})

if(MSVC)
  SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /openmp" CACHE STRING
    "Flags used by the compiler." FORCE)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_BUILD}(r${VERSION_BUILD})")
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")

set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/lib/xml")
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

# Boost Source File
find_package(BoostSrc REQUIRED)

# MPI
find_package(MPI)
IF(MPI_FOUND AND NOT ALPS_FOR_VISTRAILS)
  LINK_LIBRARIES(${MPI_LIBRARIES})
  include_directories(${MPI_INCLUDE_PATH})
  set(BUILD_BOOST_MPI TRUE)
ELSE(MPI_FOUND AND NOT ALPS_FOR_VISTRAILS)
  set(BUILD_BOOST_MPI FALSE)
ENDIF(MPI_FOUND AND NOT ALPS_FOR_VISTRAILS)
#add_definitions( -DBOOST_MPI_HOMOGENEOUS )

# Lapack
if (NOT ALPS_FOR_VISTRAILS)
  find_package(Lapack)
  if(MSVC)
    file(COPY ${LAPACK_LIBRARY} ${BLAS_LIBRARY} DESTINATION ${PROJECT_BINARY_DIR}/lib)
  endif(MSVC)
endif (NOT ALPS_FOR_VISTRAILS)

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
if(WIN32 AND NOT UNIX)
  find_package(PyInstaller)
  include(UsePyInstaller)
endif(WIN32 AND NOT UNIX)
find_package(PythonInterp)
set(PYTHON_SCRIPTDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
find_package(PythonLibs)
if (ALPS_FOR_VISTRAILS)
  MESSAGE(WARNING "Changing Python libraries to those of Vistrails")
  if(WIN32 AND NOT UNIX)
    SET(PYTHON_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Python25/include)
    SET(PYTHON_LIBRARY "${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Python25/libs/python25.lib")
  ELSE(WIN32 AND NOT UNIX)
    SET(PYTHON_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Contents/Frameworks/Python.framework/Headers)
    SET(PYTHON_LIBRARY "-F${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/Contents/Frameworks -framework Python")
  endif(WIN32 AND NOT UNIX)
  set(PYTHONLIBS_FOUND TRUE)
endif(ALPS_FOR_VISTRAILS)

IF (PYTHONLIBS_FOUND)
  SET(ALPS_HAVE_PYTHON 1)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS}) 
  link_libraries(${PYTHON_LIBRARY})
ELSE (PYTHONLIBS_FOUND)
  set(BUILD_BOOST_PYTHON FALSE)
ENDIF (PYTHONLIBS_FOUND)

######################################################################
# configuration files
######################################################################
configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/src/alps/config.h)

configure_file(lib/xml/lattices.xml.in ${PROJECT_BINARY_DIR}/lib/xml/lattices.xml)
configure_file(lib/xml/models.xml.in ${PROJECT_BINARY_DIR}/lib/xml/models.xml)
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
# install dependent libraries under Windows, from applications
# we'll need something similar for MacOS X
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
            DESTINATION lib
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
    install(FILES ${PROJECT_BINARY_DIR}/lib/${SQLite_NAME}.lib DESTINATION lib)
    install(FILES ${PROJECT_BINARY_DIR}/bin/${SQLite_NAME}.dll DESTINATION bin)
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
    install(FILES ${PROJECT_BINARY_DIR}/lib/${LPSolve_NAME}.lib DESTINATION lib)
    install(FILES ${PROJECT_BINARY_DIR}/bin/${LPSolve_NAME}.dll DESTINATION bin)
  endif(LPSolve_LIBRARIES)
  # alps & boost
  add_custom_target(libalps ALL
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
    COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib/xml
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/archive.exe ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/parameter2xml.exe ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/ising.exe ${PROJECT_BINARY_DIR}/bin
    COMMAND cmake -E copy_directory ${CMAKE_INSTALL_PREFIX}/lib/xml ${PROJECT_BINARY_DIR}/lib/xml
  )
endif(WIN32 AND NOT UNIX AND NOT ALPS_FOR_VISTRAILS)



######################################################################
# installation
######################################################################
include(InstallRequiredSystemLibraries)

if (ALPS_INSTALL_HEADERS)
  set(ALPS_HEADER_DIR "include")
  install(DIRECTORY src/alps src/boost DESTINATION ${ALPS_HEADER_DIR}
          FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h DESTINATION ${ALPS_HEADER_DIR}/alps)

  if(ALPS_PACKAGE_LIBRARIES)
     install(DIRECTORY ${Boost_ROOT_DIR}/boost DESTINATION ${ALPS_HEADER_DIR}
             FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
             PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
     install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/include DESTINATION .)
  endif(ALPS_PACKAGE_LIBRARIES)

endif(ALPS_INSTALL_HEADERS)



install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION lib)
if(ALPS_PACKAGE_LIBRARIES)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/lib DESTINATION .)
endif(ALPS_PACKAGE_LIBRARIES)

if (ALPS_FOR_VISTRAILS)
  install(DIRECTORY lib/pyalps DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}
        FILES_MATCHING PATTERN "*.py"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)
  install(DIRECTORY vistrails/userpackages/alps DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PACKAGE_DIR}
        FILES_MATCHING PATTERN "*.py"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  install(DIRECTORY lib/xml DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR}
        FILES_MATCHING PATTERN "*.xsl"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
else (ALPS_FOR_VISTRAILS)
  install(DIRECTORY lib/xml DESTINATION lib
        FILES_MATCHING PATTERN "*.xsl"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION lib)

  install(FILES config/CMake/UseALPS.cmake ${PROJECT_BINARY_DIR}/config/CMake/ALPSConfig.cmake
        DESTINATION share/alps)

  install(DIRECTORY lib/pyalps DESTINATION lib
        FILES_MATCHING PATTERN "*.py"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
endif (ALPS_FOR_VISTRAILS)

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
  
##########################################################################
# Binary packages                                                        #
##########################################################################

if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_NAME "ALPS Vistrails Patch")
else(ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_NAME "ALPS")
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
set(CPACK_PACKAGING_INSTALL_PREFIX "/")
if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "VisTrails")
else (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "ALPS")
endif (ALPS_FOR_VISTRAILS)
set(CPACK_PACKAGE_DEFAULT_LOCATION ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_VERSION_PATCH "${ALPS_VERSION_BUILD}-${VERSION_BUILD}")

if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  # NOTE: No ALPS icon yet
#  set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/tools/build/CMake\\\\InstallIcon.bmp")
#  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}$")
  set(CPACK_NSIS_CONTACT "wistaria@comp-phys.org")
  set(CPACK_NSIS_HELP_LINK "http://alps.comp-phys.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http://alps.comp-phys.org")
  set(CPACK_NSIS_MODIFY_PATH ON)
  
  # Encode the compiler name in the package 
  if(NOT ALPS_FOR_VISTRAILS)
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
  endif(NOT ALPS_FOR_VISTRAILS)
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

if (NOT ALPS_FOR_VISTRAILS)
  ENABLE_TESTING()
  add_subdirectory(test)
endif (NOT ALPS_FOR_VISTRAILS)
