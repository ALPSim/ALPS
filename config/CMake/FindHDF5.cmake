#
# this module look for HDF5 (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# HDF5_INCLUDE_DIR  = where hdf5.h can be found
# HDF5_LIBRARY      = the library to link against (hdf5 etc)
# HDF5_FOUND        = set to true after finding the library
#
IF(EXISTS ${PROJECT_CMAKE}/Hdf5Config.cmake)
  INCLUDE(${PROJECT_CMAKE}/Hdf5Config.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/Hdf5Config.cmake)

IF(Hdf5_INCLUDE_DIRS)

  FIND_PATH(HDF5_INCLUDE_DIR hdf5.h ${Hdf5_INCLUDE_DIRS})
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${Hdf5_LIBRARY_DIRS})

ELSE(Hdf5_INCLUDE_DIRS)

  SET(TRIAL_LIBRARY_PATHS
    $ENV{PHDF5_HOME}/lib
    $ENV{HDF5_HOME}/lib
    $ENV{HDF_HOME}/lib
    /usr/apps/lib
    /usr/lib 
    /opt/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/hdf5/lib
    /sw/lib
    $ENV{HOMEDRIVE}/HDF5/lib
    )

  SET(TRIAL_INCLUDE_PATHS
    $ENV{PHDF5_HOME}/include
    $ENV{HDF5_HOME}/include
    $ENV{HDF_HOME}/include
    /usr/apps/include
    /usr/include
    /opt/include
    /opt/local/include
    /usr/local/include
    /usr/local/hdf5/include
    /sw/include
    $ENV{HOMEDRIVE}/HDF5/include
    )

  IF($ENV{HDF5_DIR} MATCHES "hdf")
    MESSAGE(STATUS "Using environment variable HDF5_DIR.")
    SET(TRIAL_LIBRARY_PATHS $ENV{HDF5_DIR}/lib ${TRIAL_LIBRARY_PATHS} )
    SET(TRIAL_INCLUDE_PATHS $ENV{HDF5_DIR}/include ${TRIAL_INCLUDE_PATHS} )
  ENDIF($ENV{HDF5_DIR} MATCHES "hdf")
  
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${TRIAL_LIBRARY_PATHS})
  FIND_PATH(HDF5_INCLUDE_DIR hdf5.h ${TRIAL_INCLUDE_PATHS} )

ENDIF(Hdf5_INCLUDE_DIRS)

IF(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
  SET(HDF5_FOUND 1 CACHE BOOL "Found hdf5 library")
ELSE(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
  SET(HDF5_FOUND 0 CACHE BOOL "Not fount hdf5 library")
ENDIF(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)

MARK_AS_ADVANCED(
  HDF5_INCLUDE_DIR 
  HDF5_LIBRARY 
  HDF5_FOUND
)
