# Contribution guidance

[CONTRIBUTING.md](CONTRIBUTING.md) is the contributor guide (build, process,
provenance and scientific credit), [tests/README.md](tests/README.md) the
testing guide, and [REVIEW.md](REVIEW.md) what reviewers look for. The rules
below are the ones AI-assisted changes most often miss.

- Anything a person will read must account for limited human attention.
- Build and run the tests that cover the change before opening a pull request
  ([run the tests](CONTRIBUTING.md#run-the-tests)). In the pull request's
  "Tested" section, name the commit you built and quote the command and what it
  reported, not a claim that it passes.
- Add only the test the change needs
  ([when to add a test](tests/README.md#when-to-add-a-test)). A bug fix adds
  the one test that fails without the fix, and the pull request says so.
- For a change that affects numerical results, fill in "Validated against"
  with the comparison you made, not a description of one.
- Record the source of copied, translated or adapted code next to it, keep
  copyright and license notices, and flag uncertain provenance for maintainer
  review. Do not invent attribution or claim independence without evidence.
- For scientific credit, update `CITATION.cff` and `CITATIONS.yaml` per
  [the citation instructions](.github/scripts/citations/README.md).
  `CITATION.md` is generated; do not edit it.
