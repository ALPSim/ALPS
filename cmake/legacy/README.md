# Legacy CMake

ALPS 2.x consumer interfaces, kept as deprecated layers over the exported SDK
targets. `ALPSLegacyInterfaces.cmake` installs them next to `ALPSConfig.cmake`:

- `UseALPS.cmake`: `include(${ALPS_USE_FILE})`, `${ALPS_LIBRARIES}`, `add_alps_test`
- `UsePyALPS.cmake`: `alps_target_link_pyalps`
- `include.mk.in`: the Makefile fragment `share/alps/include.mk`

`tests/cmake/legacy/` builds a CMake and a Makefile consumer against them.

To remove them, delete this directory and `tests/cmake/legacy/`, the
`include(cmake/legacy/ALPSLegacyInterfaces.cmake)` line in the top-level
`CMakeLists.txt`, the `ALPS_USE_FILE` and `ALPS_PYTHON_USE_FILE` lines in
`cmake/ALPSConfig.cmake.in`, and the `test_legacy_*` tests in
`tests/cmake/test_sdk_contracts.py`.

Unused ALPS 2.x find modules, test macros, Boost build helpers and installer
scripts have been removed. The consumer interfaces described above remain.
