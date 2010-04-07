######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything 
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
set(FOUND_CXXENV 0)
include(TestCXXAcceptsFlag)
include(CheckCCompilerFlag)

######################################################
# Try to identify CPU identity
######################################################
set(CPU_IDENTITY "generic")
include(CheckProcessorID)

#------------------------------------
# On Jaguar use CrayXT.cmake
#------------------------------------
IF($ENV{HOST} MATCHES "jaguar")
  MESSAGE("  Working on jaguar.")	 
  include(CrayXT)
  set(FOUND_CXXENV 1)
  set(CPU_IDENTITY "barcelona")
ENDIF($ENV{HOST} MATCHES "jaguar")

#------------------------------------
# Check if using IBM compilers
#------------------------------------
IF($ENV{CXX} MATCHES "xlC")
  include(IBMCompilers)
  set(FOUND_CXXENV 1)
ENDIF($ENV{CXX} MATCHES "xlC")

#------------------------------------
# Check if using Intel compilers
#------------------------------------
IF($ENV{CXX} MATCHES "icpc")
  include(IntelCompilers)
  set(FOUND_CXXENV 1)
ENDIF($ENV{CXX} MATCHES "icpc")

#------------------------------------
# other compilers, e.g., mpicxx 
# most likely *unix with gnu or intel compilers
# using "-restrict" option to find out if intel compilers are backend.
#------------------------------------
IF(NOT FOUND_CXXENV)
  IF(CMAKE_COMPILER_IS_GNUCXX)
    include(GNUCompilers)
  ELSE(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_TRY_INTEL_CXX_FLAGS "-restrict")
    CHECK_CXX_ACCEPTS_FLAG(${CMAKE_TRY_INTEL_CXX_FLAGS} INTEL_CXX_FLAGS)
    IF(INTEL_CXX_FLAGS)
      include(IntelCompilers)
      set(FOUND_CXXENV 1)
    ENDIF(INTEL_CXX_FLAGS)
  ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ENDIF(NOT FOUND_CXXENV)
