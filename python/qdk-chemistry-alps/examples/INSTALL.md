# Reproduce the QDK–ALPS examples

This development snapshot includes the complementary-operator molecular solver,
classic time evolution, the QDK adapter, and the advanced examples. It does not
require access to the private source repository named in the provenance comments.

The instructions below target **macOS Apple Silicon with Python 3.13**. The
ALPS SDK, its Python bindings, and the pinned QDK core are built from source.
Preinstalled third-party native dependencies can be reused. Linux/WSL users
should use the pinned QDK checkout's installation guide for platform-specific
native dependencies; the local verification does not establish Linux support.

## 1. Get the public ALPS branch and create an isolated environment

```bash
git clone --single-branch --branch examples/qdk-complementary-dmrg https://github.com/ALPSim/ALPS.git
cd ALPS
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
unset PYTHONPATH QDK_CHEMISTRY_RESOURCES_PATH
```

Record `git rev-parse HEAD` with results. Do not import an older checkout or
copy Python files over an installed QDK binary. All following commands use
this same virtual environment.

## 2. Build the ALPS SDK and bindings

Install Xcode command-line tools if needed. Prebuilt Homebrew dependencies
avoid rebuilding Boost, HDF5 and Eigen:

```bash
brew install cmake ninja boost hdf5 eigen openblas
export CMAKE_PREFIX_PATH="$(brew --prefix boost):$(brew --prefix hdf5):$(brew --prefix eigen):$(brew --prefix openblas)"
cmake --preset sdk
cmake --build --preset sdk --target install --parallel 2
export ALPS_DIR="$PWD/_build/sdk/install/share/alps"
PYALPS_BUNDLE_APPLICATIONS=OFF CMAKE_BUILD_PARALLEL_LEVEL=2 \
  python -m pip install ./python/pyalps \
    -c python/qdk-chemistry-alps/examples/requirements-reproduce.txt \
    -C cmake.define.PYALPS_BUILD_SOLVERS=OFF \
    -C build-dir="$PWD/_build/python"
```

The SDK preset builds the models, exact-diagonalization and classic-DMRG
libraries needed here. The optional application executables and unrelated
Monte Carlo solver bindings are unnecessary for these examples. Keep the
build directories for incremental rebuilds.

## 3. Install the compatible QDK MPS snapshot

Use QDK Chemistry PR [#540](https://github.com/microsoft/qdk-chemistry/pull/540)
at `fb286eb2d9dde51d960f6953cb726cb45d52c12b`. It already includes the container
changes from PR #751 at `079029fad1b706a0d59aa3f1ae603aa26e2ba913`.
There is only one QDK checkout and no cherry-pick or Python overlay.

From the ALPS repository root:

```bash
git init ../qdk-chemistry
git -C ../qdk-chemistry remote add origin https://github.com/microsoft/qdk-chemistry.git
git -C ../qdk-chemistry fetch --depth 1 origin fb286eb2d9dde51d960f6953cb726cb45d52c12b
git -C ../qdk-chemistry checkout --detach FETCH_HEAD
cd ../qdk-chemistry
```

Set `QDK_DEPS` to a compatible existing QDK **third-party dependency prefix**.
Microsoft developers can reuse the prefix prepared for this checkout. The
source build below explicitly excludes any old QDK library in that prefix.
If native dependencies are not installed, the pinned repository provides its
own manifest-driven CI helper:

```bash
brew install nlohmann-json wget
export QDK_DEPS="$PWD/../qdk-deps"
INSTALL_PREFIX="$QDK_DEPS" BUILD_DIR="$PWD/../qdk-deps-build" \
  MARCH=armv8-a JOBS=2 \
  bash .pipelines/install-scripts/install-cpp-deps.sh \
    "$PWD/cpp/manifest/qdk-chemistry/cgmanifest.json" \
    "$PWD/external/macis/manifest/cgmanifest.json" auto
```

That optional helper builds substantial dependencies, including Libint, and
can take considerably longer than the examples. It was inspected against
upstream CI, not rerun as part of the local verification. Reuse an existing
compatible dependency installation where available. The validation record
separately identifies the versions actually used.

Compile and install the pinned QDK core and Python bindings:

```bash
CMAKE_PREFIX_PATH="$QDK_DEPS:$(brew --prefix)" \
CMAKE_BUILD_PARALLEL_LEVEL=2 QDK_ENABLE_OPENMP=OFF QDK_UARCH=armv8-a \
  python -m pip install ./python \
    -c ../ALPS/python/qdk-chemistry-alps/examples/requirements-reproduce.txt \
    -C build-dir=build/python \
    -C cmake.define.CMAKE_DISABLE_FIND_PACKAGE_qdk=TRUE \
    -C cmake.define.QDK_ALLOW_DEPENDENCY_FETCH=OFF \
    -C cmake.define.CMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF \
    -C cmake.define.BUILD_TESTING=OFF \
    -C cmake.define.MACIS_ENABLE_TESTS=OFF \
    -C cmake.define.CMAKE_CXX_FLAGS="-I$QDK_DEPS/include"
python -c 'from qdk_chemistry.data import MPSContainer, MPSSite; import qdk_chemistry; print(qdk_chemistry.__version__)'
cd ../ALPS
```

The include flag keeps headers consistent with the prefix's spdlog library.
The pinned Python metadata selects QDK runtime `1.32.3` and NumPy `>=2,<2.4`.
The checked-in constraints file pins the scientific Python versions used in
the clean verification. Respect the QDK constraints. `RdmContainer` is a separate optional QDK capability;
the MPS examples require `MPSContainer` and `MPSSite` instead.

## 4. Install the adapter and run the examples

```bash
python -m pip install './python/qdk-chemistry-alps[examples,test]' \
  -c python/qdk-chemistry-alps/examples/requirements-reproduce.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QSHARP_PYTHON_TELEMETRY=false
python -m pytest python/qdk-chemistry-alps --require-solver-stack -q
cd python/qdk-chemistry-alps/examples
python quickstart.py
python spin_quench.py --output results/spin
python nitrogen/run.py prepare --distance 2.5 --folder results/n2
python nitrogen/run.py solve --folder results/n2 --bonds 8 16 32
python render.py results/n2/solve.json
```

The six-site Hubbard quickstart compares exact diagonalization and DMRG.
The spin quench compares native time evolution, Q# execution and an independent
exact reference. Nitrogen prepares a full CAS(6e,6o) Hamiltonian and checks the
final MPS, signed coefficients, energy/RDM consistency and serialization.
The molecular example is small enough for the independent exact reference.

Optional preparation and resource estimation:

```bash
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py circuit --folder results/n2 --bond 16 --bits 8
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py qpe_resources --folder results/n2 --bond 16 --bits 8
```

These have different meanings: `circuit` simulates preparation before reset,
including all bond outcomes; `qpe_resources` counts a coarse QPE demonstration
and does not execute QPE. Higher precision may require much more time or
memory. See [nitrogen/README.md](nitrogen/README.md) and the
[example index](README.md) for the remaining supported demonstrations.

## Verification scope

See [validation/README.md](validation/README.md). The clean-install check uses
fresh ALPS/QDK source checkouts and a fresh virtual environment, with rebuilt
native ALPS and QDK code and preinstalled third-party native dependencies.
It is not a claim that a blank machine's dependency provisioning was retested.
Earlier figure measurements remain labeled with their original runtime.
