# CI change detection

The PR workflow checks six areas defined in `areas.json`. `fingerprint.py`
compares the complete PR diff against its base (master pushes use the previous
commit), then hashes each area's inputs, configuration, and runner image.
Successful jobs save an exact-key pass marker; compiler caches are separate
and never establish that tests passed.

Comment-only changes can reuse native or packaging evidence. Static checks
compare raw bytes. Unsupported or ambiguous syntax falls back to raw bytes;
unknown paths, file additions/deletions, missing history, and Git errors cause
checks to run. Citation authorities and generated snapshots are raw global
inputs. The manual `force` input bypasses pass markers and area selection.

Run the detector's tests without external dependencies:

```sh
python3 -m unittest ci/test_fingerprint.py
python3 ci/fingerprint.py --check-workflow .github/workflows/ci.yml
```

The `mpi` CMake preset registers retained MPI executables with CTest through
`register_mpi_tests.cmake`. Fixtures retain their original stdin/golden-output
files. Tests use separate working directories and the rank counts required
by their algorithms; the process-group fixture requires eight ranks.

This pipeline is adapted from [skilledwolf's PR #166 branch](https://github.com/skilledwolf/ALPS/tree/b5e6f650408bafd8fba77a3dfcc04426f009f4fa).
Test paths, MPI registration and sanitizer flags are adapted to the SDK branch's
retained native test framework. Compatibility runs weekly, manually and before
release publication; see [CI coverage](../CONTRIBUTING.md#ci-coverage).
