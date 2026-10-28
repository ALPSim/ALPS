# Contribution guidance

Follow [CONTRIBUTING.md](CONTRIBUTING.md) for contribution requirements,
including its provenance and scientific credit policy, and
[REVIEW.md](REVIEW.md) for what reviewers look for and what a test is for.

## Verify before opening a pull request

Build and run the tests that cover the change, and report the actual output
rather than a claim that it passes. `CMakePresets.json` defines the presets;
`dev` builds the C++ library only, `default` adds the applications, `mpi`
enables MPI.

```sh
cmake --preset dev
cmake --build --preset dev --parallel 2
ctest --preset dev -R <component or test pattern>
```

Application and tutorial tests need the `default` preset. Python changes are
tested after the editable install described in
[CONTRIBUTING.md](CONTRIBUTING.md#build):

```sh
python -m pytest python/pyalps/tests -q
```

For a change that affects numerical results, compare against an exact or
independent value (exact diagonalization, an analytic result, a previous
release) and put that comparison in the pull request's "Validated against"
section. Add only tests that would fail for a real bug; do not add tests that
pin internal details or API shapes that may change.

## Provenance and credit

When copying, translating, or substantially adapting external code, record
the source and nature of the adaptation near the affected code. Preserve
copyright and license notices, and flag uncertain provenance or licensing
for maintainer review. Do not invent attribution or claim independence
without evidence.

For scientific credit changes, update `CITATION.cff` and `CITATIONS.yaml`
using [the citation maintenance instructions](.github/scripts/citations/README.md).
`CITATION.md` is generated; do not edit it directly.
