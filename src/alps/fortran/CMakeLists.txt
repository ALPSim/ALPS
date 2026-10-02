# This is a C++ bridge; only its consumers need a Fortran compiler.
add_library(alps_fortran STATIC src/fwrapper_impl.C)
add_library(ALPS::fortran ALIAS alps_fortran)
set_target_properties(alps_fortran PROPERTIES EXPORT_NAME fortran)
target_link_libraries(alps_fortran PUBLIC alps)
# The F77 entry points accept several argument types.
target_compile_options(alps_fortran INTERFACE
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch>")

target_sources(alps_headers INTERFACE
  FILE_SET fortran_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/fortran/alps_fortran.h
    include/alps/fortran/fortran_wrapper.h
    include/alps/fortran/fwrapper_impl.h)

install(TARGETS alps_fortran EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

alps_register_module(NAME fortran HEADER_SET fortran_headers
  TARGETS alps alps_fortran
  DEPENDS legacy_parameters parapack
  PRIVATE_FILES
    src/fwrapper_impl.C
)
