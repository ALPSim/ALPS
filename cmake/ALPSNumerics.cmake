# SPDX-License-Identifier: MIT
# Preserve provider targets for the default build; delegate explicit provider,
# linkage and integer-ABI requests to CMake's standard numerical finders.
set(ALPS_BLA_VENDOR "${BLA_VENDOR}")
set(ALPS_BLA_STATIC "${BLA_STATIC}")
if(BLA_F95)
  message(FATAL_ERROR "ALPS uses the BLAS/LAPACK F77 ABI; BLA_F95 is not supported")
endif()
if(NOT DEFINED BLA_SIZEOF_INTEGER)
  set(BLA_SIZEOF_INTEGER 4)
endif()
if(NOT BLA_SIZEOF_INTEGER MATCHES "^(4|8)$")
  message(FATAL_ERROR "ALPS requires BLA_SIZEOF_INTEGER=4 or 8 (not ANY)")
endif()
set(BIND_FORTRAN_INTEGER_8 OFF)
if(BLA_SIZEOF_INTEGER EQUAL 8)
  set(BIND_FORTRAN_INTEGER_8 ON)
endif()

if(BLA_SIZEOF_INTEGER EQUAL 4
   AND NOT BLA_STATIC
   AND NOT BLA_PREFER_PKGCONFIG
   AND (NOT BLA_VENDOR OR BLA_VENDOR STREQUAL "All"))
  find_package(OpenBLAS CONFIG QUIET)
  if(OpenBLAS_FOUND AND TARGET OpenBLAS::OpenBLAS)
    # Config packages do not consume BLA_SIZEOF_INTEGER. Check their public ABI.
    include(CheckCXXSourceCompiles)
    include(CMakePushCheckState)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_LIBRARIES OpenBLAS::OpenBLAS)
    unset(ALPS_OPENBLAS_LP64 CACHE)
    check_cxx_source_compiles(
      "#include <openblas_config.h>
      static_assert(sizeof(blasint) == 4, \"OpenBLAS must use LP64 integers\");
      int main() {}" ALPS_OPENBLAS_LP64)
    cmake_pop_check_state()
    if(NOT ALPS_OPENBLAS_LP64)
      message(FATAL_ERROR "The OpenBLAS config package does not provide the requested LP64 ABI")
    endif()
    set(ALPS_BLAS_TARGET OpenBLAS::OpenBLAS)
    set(ALPS_BLAS_PROVIDER OpenBLAS)
  endif()
  # NAMES keeps config mode even with vcpkg's LAPACK wrapper.
  find_package(lapack CONFIG QUIET NAMES lapack)
  if(lapack_FOUND AND TARGET lapack)
    set(ALPS_LAPACK_TARGET lapack)
    set(ALPS_LAPACK_PROVIDER lapack)
  endif()
endif()

# Numerical templates and the application suite share one explicit ABI.
if(NOT ALPS_BLAS_TARGET)
  find_package(BLAS REQUIRED)
  set(ALPS_BLAS_TARGET BLAS::BLAS)
  set(ALPS_BLAS_PROVIDER BLAS)
endif()
if(NOT ALPS_LAPACK_TARGET)
  find_package(LAPACK REQUIRED)
  set(ALPS_LAPACK_TARGET LAPACK::LAPACK)
  set(ALPS_LAPACK_PROVIDER LAPACK)
endif()
# Some dependency wrappers modify these variables while finding transitive
# libraries. Keep the user's selection stable for the installed SDK.
if(ALPS_BLA_VENDOR
   AND NOT ALPS_BLA_VENDOR STREQUAL "All"
   AND NOT BLA_VENDOR STREQUAL ALPS_BLA_VENDOR)
  message(
    FATAL_ERROR
      "A dependency finder overrode BLA_VENDOR=${ALPS_BLA_VENDOR}; use a toolchain that honors the requested numerical provider"
  )
endif()
if(ALPS_BLA_STATIC AND NOT BLA_STATIC)
  message(
    FATAL_ERROR
      "A dependency finder overrode BLA_STATIC=ON; use a toolchain with static numerical libraries")
endif()
set(BLA_VENDOR "${ALPS_BLA_VENDOR}")
set(BLA_STATIC "${ALPS_BLA_STATIC}")
set(ALPS_HAVE_BLAS 1)
set(ALPS_HAVE_LAPACK 1)
