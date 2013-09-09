# Source: http://ericscottbarr.com/blog/2012/03/sphinx-and-cmake-beautiful-documentation-for-c-projects
# Adapted for ALPS by
#                      Michele Dolfi <dolfim@phys.ethz.ch>

find_package(Sphinx REQUIRED)
find_package(Doxygen REQUIRED)

## Doxygen settings
# location of doxygen xml output
set(DOXYGENDOCS_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

#configure_file(get_gasp.in "${CMAKE_CURRENT_BINARY_DIR}/get_gasp")
#
#add_custom_target(gasp ALL
#    "/bin/bash"
#    "${CMAKE_CURRENT_BINARY_DIR}/get_gasp"
#    COMMENT "Getting Gasp")

configure_file(doxygen.conf.in "${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf")

## Sphinx settings
# path to ext modules
set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# configured documentation tools and intermediate build results
set(SPHINX_BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

configure_file(conf.py.in "${CMAKE_CURRENT_BINARY_DIR}/conf.py")

## Doxygen target
add_custom_target(doxygendocs ALL
    ${DOXYGEN_EXECUTABLE}
    "${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf"
    COMMENT "Building Doxygen XML documentation for Sphinx")


set(PYALPS_DEPS pyalea_c pymcdata_c pytools_c  pyngsparams_c pyngshdf5_c pyngsbase_c  
    pyngsobservable_c pyngsobservables_c pyngsresult_c pyngsresults_c pyngsapi_c pyngsrandom01_c
)

if(LAPACK_FOUND)
   set(PYALPS_DEPS ${PYALPS_DEPS} maxent_c dwa_c)
endif(LAPACK_FOUND)

message(STATUS "Sphinx dependencies: ${PYALPS_DEPS}")
## Sphinx targets
add_custom_target(html ALL
    ${SPHINX_EXECUTABLE}
    -b html
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    -d "${SPHINX_CACHE_DIR}"
    "${SPHINX_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
    COMMENT "Building HTML documentation with Sphinx")
add_dependencies(html doxygendocs gasp ${PYALPS_DEPS})

add_custom_target(man ALL
    ${SPHINX_EXECUTABLE}
    -b man
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    -d "${SPHINX_CACHE_DIR}"
    "${SPHINX_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/man"
    COMMENT "Building manual pages with Sphinx")
add_dependencies(man doxygendocs gasp)



