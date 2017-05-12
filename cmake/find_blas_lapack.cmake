# Search for BLAS and LAPACK unless told otherwise

if (BLA_VENDOR MATCHES "DEFAULT")
   set(BLAS_LIBRARIES)
   set(LAPACK_LIBRARIES)
   message(STATUS "Search for BLAS and LAPACK skipped.")
elseif (BLAS_USE_MKL)
   if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
      set(BLAS_LIBRARIES)
      set(LAPACK_LIBRARIES)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")
      message(STATUS "Search for BLAS and LAPACK skipped, using MKL instead")
   else()
      message(FATAL_ERROR 
      "Do not know how to use MKL with non-Intel compiler."
      " Try setting `-DBLA_VENDOR=Intel` instead (see `cmake --help-module=FindBLAS`).")
   endif()
else()
   find_package(BLAS REQUIRED)
   find_package(LAPACK REQUIRED)
endif()
