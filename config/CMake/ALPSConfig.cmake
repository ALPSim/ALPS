#-----------------------------------------------------------------------------
#
# ALPSConfig.cmake - ALPS CMake configuration file for external projects.
#
# This file is configured by ALPS and used by the UseALPS.cmake module
# to load ALPS's settings for an external project.

# The ALPS include file directories.
set(ALPS_INCLUDE_DIRS "/Users/gamperl/ALPS/include;/Users/gamperl/boost_1_40_0")

# The ALPS library directories.
set(ALPS_LIBRARY_DIRS "/Users/gamperl/ALPS/lib")

# The ALPS runtime library directories.  Note that if
# ALPS_CONFIGURATION_TYPES is set (see below) then these directories
# will be the parent directories under which there will be a directory
# of runtime binaries for each configuration type.
set(ALPS_RUNTIME_LIBRARY_DIRS "/Users/gamperl/ALPS/bin")

# The C and C++ flags added by ALPS to the cmake-configured flags.
set(ALPS_REQUIRED_C_FLAGS "")
set(ALPS_REQUIRED_CXX_FLAGS "")
set(ALPS_REQUIRED_EXE_LINKER_FLAGS "")
set(ALPS_REQUIRED_SHARED_LINKER_FLAGS "")
set(ALPS_REQUIRED_MODULE_LINKER_FLAGS "")

# The ALPS version number.
SET(ALPS_VERSION_MAJOR "1")
SET(ALPS_VERSION_MINOR "3")
SET(ALPS_VERSION_BUILD "3")

# The location of the UseALPS.cmake file.
set(ALPS_USE_FILE "/Users/gamperl/ALPS/share/alps/UseALPS.cmake")

# ALPS Configuration options.
set(ALPS_BUILD_SHARED_LIBS "ON")
set(ALPS_BUILD_TYPE Release)
set(ALPS_CONFIGURATION_TYPES )

# MPI
set(MPI_FOUND "TRUE")
set(MPI_COMPILE_FLAGS "-D_REENTRANT")
set(MPI_EXTRA_LIBRARY "/usr/lib/libmpi.dylib;/usr/lib/libopen-rte.dylib;/usr/lib/libopen-pal.dylib")
set(MPI_INCLUDE_PATH "/usr/include")
set(MPI_LINK_FLAGS "-Wl,-u,_munmap -Wl,-multiply_defined,suppress")

# LAPACK
set(LAPACK_FOUND "")
set(LAPACK_DEFINITIONS "")
set(LAPACK_LINKER_FLAGS "")
set(LAPACK_LIBRARIES "")

# SQLite
set(SQLite_FOUND "")
set(SQLite_INCLUDE_DIR "")
set(SQLite_LIBRARIES "")
set(SQLite_DLLS "")

# LPSolve
set(LPSolve_FOUND "0")
set(LPSolve_INCLUDE_DIR "LPSolve_INCLUDE_DIR-NOTFOUND")
set(LPSolve_LIBRARIES "LPSolve_LIBRARIES-NOTFOUND")
set(LPSolve_DLLS "")
