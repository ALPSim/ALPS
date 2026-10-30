# Reviewing changes to ALPS

This file is for anyone reviewing a pull request, human or AI-assisted. It says
what matters in this repository, so that review time goes to the failures that
cost the most. [CONTRIBUTING.md](CONTRIBUTING.md) describes the process;
[tests/README.md](tests/README.md) describes the test conventions.

## What failure looks like here

ALPS output ends up in papers. A green CI run does not mean the numbers are
right. Look first for:

- Results that are wrong while the tests pass: a dropped term, a wrong sign or
  prefactor, a boundary condition that is off by one site, a sum that stops a
  step early.
- Random numbers drawn once where they should be drawn per sample, per bond or
  per site, for example when a parameter expression such as `random()` is
  evaluated at load time.
- Input that is accepted and then silently ignored: unknown parameter keys,
  trailing blocks in a parameter file, unsupported model terms that fall
  through without an error.
- Checkpoints and restarts that resume to a different state than an
  uninterrupted run would reach.
- Changes to on-disk formats (HDF5 layouts, XML schemas, parameter files) that
  older pyalps or existing result files cannot read, without a documented
  migration.
- MPI code paths that only work on rank 0, or that stop a run on one rank
  while the others wait.
- Provenance gaps: adapted code without a source note, missing citations for a
  new method ([CONTRIBUTING.md](CONTRIBUTING.md#provenance-and-scientific-credit)).

## Writing findings

- Tag each finding **Must fix**, **Should address**, **Nit** or **Pre-existing**.
  A pre-existing bug the pull request did not introduce goes into an issue, not
  into the pull request, unless the change makes it worse or easier to fix now.
- Name the file and line. For a wrong result, give a reproducer: the parameter
  file, the expected value and where it comes from, and the observed value.
- Keep nits to a few per review. Formatting is not a finding; ALPS has no
  single formatter and asks only that new code match its surroundings.
- Start the summary with the number of must-fix findings, and say plainly when
  there are none.
- Raise doubts about the approach early, in the first round, not as blocking
  comments on a finished pull request.

## Tests

- A review miss becomes a test. When a reviewer catches a wrong result, the fix
  adds the test that would have caught it, in the same pull request where that
  is practical.
- Tests check observable behaviour: a physical result against an exact or
  independent value, a file that round-trips, an input that is rejected. They
  do not pin how a result is computed, private data layouts, or API shapes that
  are expected to change. One test per behaviour, not one per function.
- "Too many tests" and "this test fixes an internal detail" are valid findings,
  and are common in AI-assisted changes. Ask for the tests to be reduced to the
  ones that would fail for a real bug.
- Follow [tests/README.md](tests/README.md#choosing-assertions-and-reference-data)
  for tolerances, seeds and fixtures. A tolerance widened to make a new
  platform pass needs a numerical reason in the pull request.

## What to skip

`third_party/`, the generated `CITATION.md`, and `legacy/` unless the pull
request changes it.
