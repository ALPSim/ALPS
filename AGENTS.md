# Contribution guidance

[CONTRIBUTING.md](CONTRIBUTING.md) is the contributor guide (build, process,
provenance and scientific credit), [tests/README.md](tests/README.md) the
testing guide, and [REVIEW.md](REVIEW.md) what reviewers look for.

- Anything a person will read must account for limited human attention. Judge
  what is relevant for the intended reader.
- Run the tests that cover the change before opening a pull request.
- Before opening a pull request, have the diff reviewed in a fresh agent
  context against REVIEW.md (`/code-review` where available), and fix, answer
  or dismiss its findings.
- Add only the tests a real bug would fail
  ([when to add a test](tests/README.md#when-to-add-a-test)).
- A change to numerical results is checked against an independent value.
- Code that cannot do what the input asks stops with an error. It does not
  print a warning and compute something else.
- Describe a bug without assigning blame. Saying which change introduced it is
  useful context.
- Record the source of copied, translated or adapted code next to it, keep
  copyright and license notices, and flag uncertain provenance for maintainer
  review. Do not invent attribution or claim independence without evidence.
- For scientific credit, update `CITATION.cff` and `CITATIONS.yaml` per
  [the citation instructions](.github/scripts/citations/README.md).
  `CITATION.md` is generated; do not edit it.
