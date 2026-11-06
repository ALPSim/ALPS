# Contribution guidance

[CONTRIBUTING.md](CONTRIBUTING.md) is the contributor guide (build, process,
provenance and scientific credit), [tests/README.md](tests/README.md) the
testing guide, and [REVIEW.md](REVIEW.md) what reviewers look for.

- Anything a person will read must account for limited human attention.
- Run the tests that cover the change before opening a pull request, and
  report what they printed, not that they pass.
- Add only the tests a real bug would fail
  ([when to add a test](tests/README.md#when-to-add-a-test)).
- A change to numerical results is checked against an independent value.
- Record the source of copied, translated or adapted code next to it, keep
  copyright and license notices, and flag uncertain provenance for maintainer
  review. Do not invent attribution or claim independence without evidence.
- For scientific credit, update `CITATION.cff` and `CITATIONS.yaml` per
  [the citation instructions](.github/scripts/citations/README.md).
  `CITATION.md` is generated; do not edit it.
