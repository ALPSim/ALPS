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
- The second code path. Most wrong results found in review in autumn 2026 were
  in the less common form of a feature, while the common form was right and
  tested: a custom measurement next to a built-in one, a multi-term operator
  next to a single-term one, a scalar parameter next to a vector one, an
  unconstrained state next to a constrained one, a value given on the command
  line next to one in the parameter file, a positional argument next to an
  option. When a change touches one form, run the other.
- Random numbers drawn once where they should be drawn per sample, per bond or
  per site, for example when a parameter expression such as `random()` is
  evaluated at load time.
- Input that is accepted and then silently ignored: unknown parameter keys,
  trailing blocks in a parameter file, unsupported model terms that fall
  through without an error, a file that only partly parses.
- The same name with a different meaning in another ALPS application, for
  example an average that one code normalises per site and another does not.
  A change that touches a shared parameter or observable name says which
  convention it follows.
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

- Say which commit you reviewed and what you built and ran. Pull requests
  change under review, and a number quoted without its commit cannot be
  checked later.
- Tag each finding **Blocking**, **Suggestion**, **Nit** or **Outside the
  diff**. A pre-existing bug the pull request did not introduce goes into an
  issue, not into the pull request, unless the change makes it worse or easier
  to fix now.
- Name the file and line. For a wrong result, give a reproducer: the parameter
  file, the expected value and where it comes from, and the observed value.
- Keep nits to a few per review. Formatting is not a finding; ALPS has no
  single formatter and asks only that new code match its surroundings.
- Start the summary with the number of blocking findings, and say plainly when
  there are none.
- A wrong approach is Blocking when you can name what it breaks, a Suggestion
  when you would only have built it differently. Name the alternative.

## Tests

- A review miss becomes a test. When a reviewer catches a wrong result, the fix
  adds the test that would have caught it, in the same pull request where that
  is practical.
- A test that still passes with the bug present is a blocking finding. For a
  bug fix, the pull request says that the new test fails without the fix.
- A test that fails the bar in [when to add a test](tests/README.md#when-to-add-a-test)
  is a finding, not a nit; ask for it to be removed before merge. This pattern
  is common in AI-assisted changes. For each test that stays, the author can
  name the bug it would catch.
- A tolerance widened to make a new platform pass needs a numerical reason in
  the pull request ([reference data](tests/README.md#choosing-assertions-and-reference-data)).

## What to skip

`third_party/`, the generated `CITATION.md`, and `legacy/` unless the pull
request changes it.
