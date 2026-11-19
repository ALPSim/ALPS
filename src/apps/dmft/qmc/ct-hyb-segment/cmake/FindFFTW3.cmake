# CMakeFind module to find (typically, user-installed) FFTW3.

# Input variables:
#   ENV{FFTW_ROOT} or FFTW3_ROOT pointing to directory containing FFTW3 installation
#   FFT3_INCLUDES or ENV{FFTW_INC} : (optional) guess to look for includes
#   FFTW3_LIBRARIES : (optional) location of FFTW3 libraries, if already known
#   ENV{FFTW_LINK} : (optional) directory to look for FFTW3 libraries

# Output variables:
#   FFTW3_FOUND : true if FFTW3 is found
#   FFTW3_LIBRARIES : FFTW3 libraries
#   FFTW3_INCLUDE_DIRS : Directory with FFTW3 headers
#
  
include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)

# FIXME: try to use pkg-config hints.

find_path(FFTW3_INCLUDE_DIR fftw3.h 
          PATHS ${FFTW3_ROOT} ${NFFT3_INCLUDES} 
          HINTS ENV FFTW_ROOT 
	  PATH_SUFFIXES include)
set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
if (NOT FFTW3_LIBRARIES) 
  find_library(FFTW3_LIBRARY_RELEASE fftw3 
               PATHS ${FFTW3_ROOT} 
               HINTS ENV FFTW_ROOT
               PATH_SUFFIXES lib)
  select_library_configurations(FFTW3)
endif()

find_package_handle_standard_args(FFTW3 DEFAULT_MSG 
                                  FFTW3_LIBRARIES
                                  FFTW3_INCLUDE_DIR)
