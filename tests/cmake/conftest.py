"""Dependency providers for configure-only contracts (no linkable binaries)."""
import pytest


@pytest.fixture
def hdf5_provider(tmp_path, request):
    # Dependencies extracted inside the build tree are omitted from CMake's
    # automatic install RPATH. Use separate directories to expose flattening.
    provider = tmp_path / "build" / "hdf5"
    (provider / "include").mkdir(parents=True)
    (provider / "hdf5-config-version.cmake").write_text('''
set(PACKAGE_VERSION @VERSION@)
if(NOT PACKAGE_VERSION VERSION_LESS PACKAGE_FIND_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
endif()
'''.replace('@VERSION@', request.param))
    # FindHDF5 falls back to header/library discovery if the config is too old.
    (provider / "include/hdf5.h").touch()
    (provider / "include/H5pubconf.h").write_text(
        f'#define H5_VERSION "{request.param}"\n')
    (provider / "lib").mkdir()
    for name in ("libhdf5.so", "libhdf5.dylib", "libhdf5.a", "hdf5.lib"):
        (provider / "lib" / name).touch()
    (provider / "hdf5-config.cmake").write_text('''
set(HDF5_VERSION @VERSION@)
set(HDF5_ENABLE_PARALLEL OFF)
set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include")
add_library(hdf5::hdf5-shared SHARED IMPORTED)
set_target_properties(hdf5::hdf5-shared PROPERTIES
  IMPORTED_CONFIGURATIONS "DEBUG;RELEASE"
  INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")
foreach(config IN ITEMS Debug Release)
  string(TOUPPER "${config}" upper)
  set(directory "${CMAKE_CURRENT_LIST_DIR}/${config}")
  file(MAKE_DIRECTORY "${directory}")
  set(library "${directory}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5${CMAKE_SHARED_LIBRARY_SUFFIX}")
  file(WRITE "${library}" "")
  set_target_properties(hdf5::hdf5-shared PROPERTIES
    IMPORTED_LOCATION_${upper} "${library}")
  if(WIN32)
    file(WRITE "${directory}/hdf5.lib" "")
    set_target_properties(hdf5::hdf5-shared PROPERTIES
      IMPORTED_IMPLIB_${upper} "${directory}/hdf5.lib")
  endif()
endforeach()
'''.replace('@VERSION@', request.param))
    return provider
