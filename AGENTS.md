# Contribution guidance

Follow [CONTRIBUTING.md](CONTRIBUTING.md) for contribution requirements,
including its provenance and scientific credit policy, and
[REVIEW.md](REVIEW.md) for what reviewers look for and what a test is for.

## Verify before opening a pull request

Build and run the tests that cover the change, following the
[build](CONTRIBUTING.md#build) and [test](CONTRIBUTING.md#run-the-tests)
instructions in CONTRIBUTING.md and the preset guide in
[tests/README.md](tests/README.md). Report the actual output rather than a
claim that it passes.

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
