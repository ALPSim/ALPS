# pyalps

Python applications and libraries for the [ALPS project](https://alps.comp-phys.org/). Install from PyPI in a virtual environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install pyalps
```

Use a separate environment to avoid overwriting commands from a source-installed ALPS SDK, or use the bindings-only installation described below.

Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` interface. Bundled applications are serial even with this extra installed.

## Command-line tools

Wheels include commands for the bundled applications, including `spinmc`, `loop`, `worm`, `dmrg`, and `sparsediag`, plus the tutorial tools `parameter2xml`, `printgraph`, `convert2xml`, `convert2text`, `plot2text`, `plot2gp`, `plot2xmgr`, `snap2vtk`, and `maxent`. For example:

```sh
parameter2xml simulation.in
spinmc --write-xml simulation.in.in.xml
convert2text simulation.in.task1.out.xml > results.txt
```

Download inputs from the [ALPS tutorials](https://alps.comp-phys.org/tutorials/); they are not included in wheels. Set `ALPS_XML_PATH` only to override the bundled XML resources.

For older tutorials, use `python` instead of the removed `alpspython` wrapper and `pyalps.plot` instead of `plot2mpl` or `extractmpl`.

## Building from source

Requires Python 3.10+, CMake 3.22+, Ninja, a C++17 compiler, BLAS/LAPACK, and HDF5. Reuse an installed ALPS C++ SDK, or build one with the `wheel-deps` preset. From the repository root, build the SDK and Python wheel with:

```sh
cmake --preset wheel-deps
cmake --build --preset wheel-deps

python -m pip install build
ALPS_DIR="$PWD/_build/wheel-deps/install/share/alps" \
  python -m build --wheel python/pyalps
```

Install the wheel from `python/pyalps/dist` with `python -m pip install`. It bundles the SDK's applications and libraries; shell launchers and Python helpers use those executables regardless of `PATH` or `ALPS_BIN_PATH`.

For bindings only, set `ALPS_DIR` to the existing SDK's `share/alps` directory and set the environment variable `PYALPS_BUNDLE_APPLICATIONS=OFF`:

```sh
ALPS_DIR="/path/to/sdk/share/alps" PYALPS_BUNDLE_APPLICATIONS=OFF \
  python -m pip install ./python/pyalps
```

This installs no command launchers and can safely share the SDK's prefix. Keep the SDK installed and add its `bin` to `PATH` for shell use. Python helpers use `ALPS_BIN_PATH` or the SDK recorded at build time; they do not search `PATH`.

In either mode, pass a full executable path to a Python helper to select a different application, including a source-built MPI application.

## Compatibility

- Wheels target individual CPython versions (3.10–3.14), not the stable ABI. Free-threading is unsupported; importing pyalps re-enables the GIL.
- Rebuild downstream C++ extensions against the same SDK revision, nanobind internals ABI, and C++ standard library as the wheel. The nanobind version is recorded in `pyalps.pyalps_config.NANOBIND_VERSION`.
- Legacy HDF5 signed-byte datasets without type metadata are read as booleans. Use a typed reader such as h5py when such a dataset contains integers.
- `pyalps.mpi` uses mpi4py's protocol, which is incompatible with Boost.MPI serialization. Communicators cannot be used as dictionary keys.

The package version comes from `ALPS_VERSION.txt`. For local prereleases, set `ALPS_VERSION_PRERELEASE` (for example, `beta.1`); release builds derive it from the Git tag.
