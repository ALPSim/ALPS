# Ietl source ownership.
target_sources(alps_headers INTERFACE
  FILE_SET ietl_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/ietl/bandlanczos.h
    include/ietl/bicgstabl.h
    include/ietl/cg.h
    include/ietl/complex.h
    include/ietl/fmatrix.h
    include/ietl/gmres.h
    include/ietl/ietl2lapack.h
    include/ietl/interface/blas.h
    include/ietl/interface/blitz.h
    include/ietl/interface/mtl.h
    include/ietl/interface/ublas.h
    include/ietl/interface/valarray.h
    include/ietl/inverse.h
    include/ietl/iteration.h
    include/ietl/jacobi.h
    include/ietl/jd.h
    include/ietl/krylov_wrapper.h
    include/ietl/lanczos.h
    include/ietl/matrix.h
    include/ietl/power.h
    include/ietl/rayleigh.h
    include/ietl/simple_arnoldi.h
    include/ietl/tmatrix.h
    include/ietl/traits.h
    include/ietl/vectorspace.h)

alps_register_module(NAME ietl HEADER_SET ietl_headers
  TARGETS alps
  DEPENDS configuration
)
