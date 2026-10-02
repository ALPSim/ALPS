# ALPS / ALPSCore contract probes

This standalone project compiles the same source against **one installed SDK at
a time**. The runner exchanges temporary HDF5 files between fresh processes.
Never combine the two SDK include paths or runtime implementations in one binary.
The probes characterize differences; a successful comparison does not establish
drop-in compatibility.

The runner covers 20 archive payload cases, params access/conversion behavior,
full params versus dictionary checkpoint loading, and wider scalar checkpoints.
Every typed archive read checks the value. Common params self-checks verify exact
values; cross-provider differences are recorded even when loading throws or
silently changes a value. `--expect` additionally rejects any measurement change
from an earlier report. A known defect recorded in a baseline is not a desired
contract: review and update that expectation when fixing it.

## Build and run

Use the already installed ALPS SDK and its compiler/dependencies. Install only
ALPSCore utilities/HDF5/params into a separate prefix. For the recorded comparison,
the reference was an unmodified checkout at
`7146b9e1f017938a94e5dae35d88467cc5ba7969`:

```sh
# Set these to existing absolute paths. Keep this scratch tree outside the SDK.
alps_repo="$PWD"
alps_sdk="$alps_repo/_build/dev/darwin-arm64/install"
deps="$alps_repo/.pixi/envs/default"
core_source=/absolute/path/to/pinned/ALPSCore
scratch=/absolute/path/to/reconciliation-build
eigen_headers=/opt/homebrew/include/eigen3
export PATH="$deps/bin:$PATH"

cmake -S "$core_source" -B "$scratch/core-sdk" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DALPS_CXX_STD=c++17 \
  -DCMAKE_C_COMPILER="$deps/bin/cc" -DCMAKE_CXX_COMPILER="$deps/bin/c++" \
  -DCMAKE_PREFIX_PATH="$deps" -DHDF5_ROOT="$deps" \
  -DEIGEN3_INCLUDE_DIR="$eigen_headers" -DALPS_INSTALL_EIGEN=OFF \
  -DTesting=OFF -DENABLE_MPI=OFF \
  '-DALPS_MODULES_DISABLE=accumulators;mc;gf;alea' \
  -DCMAKE_INSTALL_PREFIX="$scratch/core-sdk/install"
cmake --build "$scratch/core-sdk" --parallel 3
cmake --install "$scratch/core-sdk"

cmake -S "$alps_repo/tests/reconciliation" -B "$scratch/alps-probe" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DRECONCILIATION_PROVIDER=ALPS \
  -DALPS_DIR="$alps_sdk/share/alps" \
  -DCMAKE_PREFIX_PATH="$deps" -DHDF5_ROOT="$deps" \
  -DCMAKE_C_COMPILER="$deps/bin/cc" -DCMAKE_CXX_COMPILER="$deps/bin/c++"
cmake --build "$scratch/alps-probe" --parallel 2

cmake -S "$alps_repo/tests/reconciliation" -B "$scratch/core-probe" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DRECONCILIATION_PROVIDER=ALPSCore \
  -DCMAKE_PREFIX_PATH="$scratch/core-sdk/install;$deps" \
  -DCMAKE_C_COMPILER="$deps/bin/cc" -DCMAKE_CXX_COMPILER="$deps/bin/c++"
cmake --build "$scratch/core-probe" --parallel 2

python "$alps_repo/tests/reconciliation/compare.py" \
  --alps "$scratch/alps-probe/reconciliation_probe" \
  --alpscore "$scratch/core-probe/reconciliation_probe" \
  --output "$scratch/results.json"
```

The committed `baseline-darwin-arm64.json` is a historical measurement from
before ALPS rejected unsupported native parameter checkpoint types. Its
ALPSCore-to-ALPS extended checkpoint silently loaded float, unsigned and wide
integer values as zero; current ALPS rejects that checkpoint instead. Use
`--expect` only when comparing against a baseline appropriate for the tested
revision and platform, and review intentional behavior changes separately.

Adjust paths, compiler and executable suffix for other platforms. Do not use the
macOS baseline as a universal expectation: the wide-integer case uses native
`long`, whose width is recorded; on platforms with a narrower `long` its fixture
uses 42 instead of 2^40. Integer overflow results also depend on the implementation.
Run without `--expect` to characterize another platform/revision, then review the
report. The JSON records measurements, platform and executable hashes; its sibling
`.log` retains exception diagnostics, including why checkpoint loads failed.
Record source revisions and dependency versions alongside any new baseline.

For the recorded run both SDKs used Pixi Clang 21.1.8, Boost 1.91 and HDF5 1.14.6;
Core additionally used installed Eigen 5.0.1 headers. CMake was 4.4.3, C++17,
Release, shared libraries, MPI/OpenMP disabled. No reference source was patched,
no dependency was downloaded by its build, and Core's tests were disabled.
Pin `HDF5_ROOT` as well as `CMAKE_PREFIX_PATH` when multiple providers are installed:
mixing Homebrew include roots with Pixi Boost libraries can invalidate results.

The runner removes its generated archive fixtures automatically. Preserve useful
incremental SDKs; remove only the dedicated scratch builds when finished. The
foreign SDK is optional and is not downloaded or built by the main ALPS test suite.

## Gaps

This is a small contract baseline, not a full format/ABI audit. It does not yet
exercise schema-defined defaults/provenance, scientific-container adapters,
multidimensional shape inspection, compression/family files, concurrent access,
signal/crash cleanup, leak instrumentation, Python providers, MPI, or MaxEnt
scientific results.
