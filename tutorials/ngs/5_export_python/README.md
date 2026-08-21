# Export an ALPS simulation with nanobind

This example replaces the former Boost.Python export tutorial while retaining
the public `ALPS_EXPORT_SIM_TO_PYTHON` helper. Build it against an installed
ALPS SDK and the Python environment containing pyalps and nanobind:

```sh
cmake -S . -B build -GNinja \
  -DALPS_DIR=/path/to/alps/share/alps \
  -DPython_EXECUTABLE="$(command -v python)"
cmake --build build
PYTHONPATH="$PWD/build" python main.py
```

For an old export source, replace `BOOST_PYTHON_MODULE(name) {` with
`NB_MODULE(name, m) {`; the existing
`ALPS_EXPORT_SIM_TO_PYTHON(PythonName, SimulationClass)` call remains valid.
The helper imports pyalps' owning extension modules before registering the
derived class, so ALPS parameter, archive, observable, and result types are
shared safely through nanobind's process-wide type registry. The example also
uses `alps_target_link_pyalps`, supplied by `ALPS_PYTHON_USE_FILE`, to link a
consumer to the exact `libalps`, Boost, and HDF5 copies carried by a repaired
pyalps wheel. This is required for stateful library objects such as HDF5
handles; do not replace it with a second system HDF5 linkage.
