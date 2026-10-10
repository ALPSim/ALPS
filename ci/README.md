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

The `mpi` preset registers GoogleTest executables through
`alps_add_mpi_gtest` in `cmake/ALPSTesting.cmake`. Each executable runs at two
and three ranks, with isolated working directories and per-rank reports.
See [the MPI testing guide](../tests/mpi.md) for the tested contracts and limits.

This pipeline is adapted from [skilledwolf's PR #166 branch](https://github.com/skilledwolf/ALPS/tree/b5e6f650408bafd8fba77a3dfcc04426f009f4fa).
Compatibility runs weekly, manually and before release publication; see
[CI coverage](../CONTRIBUTING.md#ci-coverage).
