# Contribution guidance

Follow [CONTRIBUTING.md](CONTRIBUTING.md) for contribution requirements,
including its provenance and scientific credit policy, and
[REVIEW.md](REVIEW.md) for what reviewers look for and what a test is for.

## Verify before opening a pull request

Build and run the tests that cover the change. From the repository root, with
the dependencies from [CONTRIBUTING.md](CONTRIBUTING.md#prerequisites):

```sh
cmake --preset default
cmake --build --preset default --parallel 2
ctest --preset default -R <pattern>       # the tests that cover the change
ctest --preset default                    # the whole native suite
python -m pytest python/pyalps/tests -q   # after the editable install in CONTRIBUTING.md
```

The `dev` presets omit the applications. [tests/README.md](tests/README.md#local-workflows)
lists the other presets and the labels that select tests. In the pull
request's "Tested" section, name the commit you built and quote the command
and what it reported, not a claim that it passes.

For a change that affects numerical results, compare against an exact or
independent value (exact diagonalization, an analytic result, a previous
release) and put that comparison in the pull request's "Validated against"
section.

## Tests

A bug fix adds one test: the one that fails without the fix. Say so in the
pull request. A new feature adds a check of its result against an exact or
independent value, and one rejection test per new error a user can hit. A
refactor that changes no behaviour adds no tests; if the behaviour it moves had
none, add one behavioural test, not a mirror of the new structure. Do not add
tests of default values, parameters read back as written, constructor members,
`--help`, warning text or timing, and do not add tests for code the change does
not touch. [REVIEW.md](REVIEW.md#tests) is the reviewer's side of this.

## Provenance and credit

When copying, translating, or substantially adapting external code, record
the source and nature of the adaptation near the affected code. Preserve
copyright and license notices, and flag uncertain provenance or licensing
for maintainer review. Do not invent attribution or claim independence
without evidence.

For scientific credit changes, update `CITATION.cff` and `CITATIONS.yaml`
using [the citation maintenance instructions](.github/scripts/citations/README.md).
`CITATION.md` is generated; do not edit it directly.
