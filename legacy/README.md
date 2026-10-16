# Legacy interfaces

Deprecated ALPS 2.x interfaces, kept so existing downstream code keeps building
while it migrates. Everything here is wired into the build by one file,
`ALPSLegacy.cmake`, behind the `ALPS_LEGACY_API` option (default `ON`).

| Path | What it keeps working | Replacement |
| --- | --- | --- |
| `include/alps/legacy/params_path_constructor.ipp`, `src/` | `alps::params p(path)` for text parameter files | `alps::params_from_file(path)` from `<alps/ngs/params_from_file.hpp>` |
| `include/alps/ngs/detail/export_sim_to_python.hpp` | the old include spelling | `<pyalps/export_simulation.hpp>` |
| `include/alps/python/save_observable_to_hdf5.hpp` | the old include spelling | write through `alps::hdf5::archive` |
| `cmake/UseALPS.cmake`, `cmake/UsePyALPS.cmake` | `include(${ALPS_USE_FILE})`, `${ALPS_LIBRARIES}`, `add_alps_test`, `alps_target_link_pyalps` | imported `ALPS::` targets; the pyalps CMake package |
| `cmake/include.mk.in` | Makefile projects using `share/alps/include.mk` | CMake with the imported targets |

Each C++ interface warns at compile time. With `ALPS_LEGACY_API=OFF` none of
these files are built or installed, and `alps::params(path)` is a compile error
rather than an HDF5 open of a text file.

## Tests

`tests/` builds the C++ interfaces in the build tree (CTest label `legacy`).
The installed-SDK consumer in `tests/cmake/consumer` also builds
`tests/params_path_constructor.cpp` as `legacy_params_contract`. In
`tests/cmake/test_sdk_contracts.py`, `test_parameter_file_constructor_deprecation`
requires its deprecation warning, and `test_legacy_*` builds a CMake and a
Makefile consumer. These tests skip when the SDK has no legacy layer.

The focused sanitizer job configures `ALPS_LEGACY_API=OFF` and exercises the
explicit text factory and rejected implicit HDF5 conversion. Other native and
installed-SDK jobs exercise the default enabled layer. Use a fresh installation
prefix when disabling the layer; CMake install does not remove previously installed files.

## Archive

`archive/` holds code that no longer builds or has no supported use. Nothing in
it is compiled or installed; it is kept for reference and possible revival.

- `archive/ngs-accumulator/`: the experimental `ALPS_NGS_USE_NEW_ALEA`
  accumulator backend and its deprecated prototype headers, tests and example.
  It did not compile against the current accumulator features, and its tests
  were never registered. `backend-wiring.patch` restores the `#ifdef` switches
  it needed in the core libraries.

## Removing the legacy layer

Delete this directory. The remaining hooks are inert without it:
`include(legacy/ALPSLegacy.cmake OPTIONAL)` in the top-level `CMakeLists.txt`,
`include(ALPSLegacyConfig.cmake OPTIONAL)` in `cmake/ALPSConfig.cmake.in`, and
the `__has_include` check in `<alps/ngs/params.hpp>`, which then leaves
`params(path)` deleted. The `test_legacy_*` tests then skip; remove them and
`tests/cmake/legacy/` at the same time.
