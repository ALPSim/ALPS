# Reviewing changes to ALPS

This file is for anyone reviewing a pull request, human or AI-assisted. It says
what matters in this repository, so that review time goes to the failures that
cost the most. [CONTRIBUTING.md](CONTRIBUTING.md) describes the process;
[tests/README.md](tests/README.md) describes the test conventions.

Read the change as someone who will maintain it: is this the right change, in
the right place, and worth the code and interfaces it adds?

## What failure looks like here

ALPS output ends up in papers. A green CI run does not mean the numbers are
right. Look first for:

- Results that are wrong while the tests pass: a dropped term, a wrong sign or
  prefactor, a boundary condition that is off by one site, a sum that stops a
  step early, a quantity computed once that should be computed per sample.
- Numbers that only the author has checked. The examples and benchmarks in
  the pull request are the author's checks. For a change to results, check
  that someone other than the author has compared it with an exact or
  independent value.
- The less common form of a feature, such as a custom measurement beside a
  built-in one or a vector parameter beside a scalar one. When a change
  touches one form, run the other.
- Input that is accepted and then silently ignored: unknown parameter keys,
  trailing blocks in a parameter file, unsupported model terms that fall
  through without an error, a file that only partly parses.
- The same name with a different meaning in another ALPS application, for
  example an average that one code normalises per site and another does not.
  A change that touches a shared parameter or observable name says which
  convention it follows.
- The same code copied elsewhere. ALPS has near-duplicate code across
  applications. When a pull request fixes a bug, search for the copies: they
  still have it.
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

- Review the pull request merged into current master: that is what will land.
- Mark each finding as a **Must fix**, a **Recommendation** or a **Nit**. A
  pre-existing bug the pull request did not introduce goes into an issue, not
  into the pull request, unless the change makes it worse or easier to fix
  now.
- Say where each finding is: an inline comment, or the file and line in the
  review text, whichever serves the author better. For a wrong result, give a
  reproducer: the parameter file, the expected value and where it comes from,
  and the observed value.
- Keep nits to a few per review. Formatting is not a finding.
- A failing CI job counts against the pull request only if the change caused
  it. If the same job fails on master, or the runner was lost, do not raise it.
- When there is nothing to raise, approve without text.
- A wrong approach is a Must fix when you can name what it breaks, a
  Recommendation when you would only have built it differently. Name the
  alternative.

## Tests

- A review miss becomes a test. When a reviewer catches a wrong result, the fix
  adds the test that would have caught it, in the same pull request where that
  is practical.
- A test that still passes with the bug present is a Must fix.
- A test that does not meet the criteria in
  [when to add a test](tests/README.md#when-to-add-a-test) should be removed
  before merge. For each test that stays, the author can say which bug it
  would catch.
- A tolerance widened to make a new platform pass needs a numerical reason in
  the pull request ([reference data](tests/README.md#choosing-assertions-and-reference-data)).

## What to skip

`third_party/`, the generated `CITATION.md`, and `legacy/` unless the pull
request changes it.
