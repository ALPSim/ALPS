# Maintaining ALPS citations

Two files own distinct information:

- `CITATION.cff` owns bibliographic records and the preferred framework citation.
  It follows the unextended Citation File Format 1.2.0 schema. Each record has
  one ordinary `identifiers` entry with `type: other`, `description: ALPS reference key`,
  and a stable `value` such as `white1992`. Other identifiers remain available.
- `CITATIONS.yaml` owns the scientific-credit request and the application-to-paper
  relationships. Its `algorithm` and `implementation` lists reference those keys.
  Every component also cites the CFF `preferred-citation`. The optional `uses`
  list includes another component's references; this credits the conventional
  scheduler when an application uses it. Empty lists explicitly mean no reference
  is specified for that tier. Missing lists are invalid.

No selection rules are encoded in CFF `scope`, `notes`, or other free-text fields.
The YAML policy is checked against `policy.schema.json`; the generator also checks
reference identifiers, duplicate DOIs, component dependencies, and cycles. Reference
order follows the explicit lists, with a paper printed once even when it serves
multiple roles. Legal license terms remain in `LICENSE.txt`.

## Editing and validation

Use an existing Python environment where practical (Python 3.9 or newer):

```sh
python -m pip install -r .github/scripts/citations/requirements.txt
python .github/scripts/generate_citations.py --regenerate
python .github/scripts/generate_citations.py --check
python -m unittest discover -s tests/packaging -p '*citations.py' -v
```

Commit the two authority files, generated `CITATION.md`, and
`.github/scripts/citations/generated/`. CI regenerates all derived files and checks that
they match. To add an application, add its component mapping and connect
its startup path to the component key; do not put bibliographic text in C++.
Use `alps::print_citations(out, component)` for a compact citation box,
`alps::print_citation_details(out, component)` for complete references, or
`alps::print_copyright(out, component)` for the combined framework banner.
Conventional scheduler applications pass the component as the fourth argument to
`alps::scheduler::start`. Parapack applications register it with
`PARAPACK_SET_CITATION_COMPONENT`.

## Building and distributing

CMake consumes the checked-in generated data without running Python or parsing
YAML. It verifies SHA-256 checksums of the authorities, generator, schemas, and
derived files; stale data stops configuration with the regeneration command.
A source change also triggers this check on the next incremental build. This
prevents a modified catalog from silently compiling yesterday's recommendations.

Generated C++ data is compiled into the utilities component (`ALPS::utilities`); execution needs no data-file lookup
or network access. Updating an installed catalog does not alter a binary.

Python, PyYAML, and jsonschema are catalog-maintenance and validation tools, not
dependencies of normal native or wheel builds. The schemas are local, so
regeneration works offline. Optional citation tests use
`-DALPS_CITATION_PYTHON=/path/to/python`; C++ builds remain available when Python
is absent. The generated directory is marked as derived for GitHub diff display,
so scientific review can focus on the authority files and readable matrix.

Native and SDK-only installations include the authorities and generated guidance
in `share/alps`; wheels include them in `pyalps/share/alps`. The standalone pyalps
source distribution uses the catalog installed with its linked C++ SDK. Native ALPS source archives must retain the two authorities, generator,
schemas, and generated directory. Website maintainers
can consume these same files or generated Markdown from a release tag; the website
must not become a separately maintained citation list.

## Initial policy and publication status

The initial catalog was transcribed from the release 3.0 paper's citation policy,
application table, bibliography, and DMFT discussion at the revision recorded in
`CITATIONS.yaml`. These repository files now own ongoing updates. The source paper
is provenance, not a dependency of the build.

- QWL credits the two Wang–Landau papers and the quantum Wang–Landau paper as
  algorithm references, with the release 3.0 paper as implementation and framework.
- The aggregate DMFT table row is split into the driver, CT-INT, CT-HYB, and
  Hirsch–Fye, using the manuscript's solver discussion. The driver asks users to
  cite their chosen solver as well and prints the solver notice for its built-in
  Hirsch–Fye and CT-INT solvers; the standalone solver programs print their own.
  Arbitrary external solvers cannot be inferred.
- Scheduler-dependent applications also credit the scheduler implementation listed
  in the paper. Parapack/looper is not assumed to be the same scheduler.
- The release 3.0 preferred citation is the arXiv preprint (arXiv:2610.03884),
  recorded with its arXiv DOI and `status: preprint`; its author list matches the
  arXiv listing. Replace it with the journal record once published.
- Emanuel Gull confirmed the DMFT split and solver references on 2026-10-08.
  CT-HYB also credits Hafermann, Werner, and Gull, CPC 184, 1280 (2013), alongside
  CPC 182, 1078 (2011), as implementation references. The shipped solver's manual
  (`tutorials/06-hybridization/hybdoc.tex`) identifies the
  2013 paper with this version of the code. Boehnke et al., PRB 84, 075145 (2011),
  remains in the manual; CT-INT cites the 2005 PRB rather than the 2004 JETP letter.
  The improved-estimator paper's existing manual-only treatment is unchanged.

Recommendations are application-level; parameter-dependent algorithm selection is
not inferred.

## CLI printing policy

Calculation startup prints one compact, decorative citation box per invocation, on the master
rank for MPI calculations. This includes looper's stdin path. Stdin calculations
write the startup notice to stderr to preserve their numerical stdout format.
Multiple tasks in one invocation share the notice. The framework paper is always
reference [1], followed by algorithm and implementation references, with each paper
listed once. Startup uses short authors and publication coordinates; unpublished
papers retain their title and status, and preprints also show their arXiv
identifier. Up to three authors are listed; longer lists use the first author
followed by "et al.". Complete titles and DOI links appear in the detailed guidance.

Set `ALPS_NO_CITATIONS=1` to suppress automatic compact citation notices in batch
runs. Only the literal value `1` disables them: unset, empty, `0`, and other values
leave notices enabled. Copyright and license notices remain visible.
`citation_text`, `citation_details`, `print_citation_details`, and explicit
`--citations` queries always provide their usual output, including when
`ALPS_NO_CITATIONS=1`.

Each application provides these standalone queries:

- `--citations`: print complete guidance with paper titles and DOI links, and exit successfully.
- `--help` / `-h`: print usage and information options, without citation guidance.
- `--license` / `-l`: print legal license terms, without citation guidance.

Queries require no input files and create no simulation outputs. An optional `--mpi`
is accepted; in MPI builds, rank zero prints the query result even when `--mpi` is
omitted. Mixing queries or adding calculation arguments is rejected. `--` ends option
recognition. Ordinary parsing failures do not print a citation notice.

New entry points should use `alps::handle_cli_information` from `ALPS::cli` before reading input.
The query helpers retain the `<alps/utility/cli.hpp>` public include. CT-HYB and
CT-INT choose their profiles in the CLI wrappers and print startup notices in the
shared SDK solver implementations, so Python calls receive the same guidance.
The helper owns MPI initialization only when needed and preserves MPI owned by its
caller. The scheduler, mcoptions, and both parapack implementations use this helper.
`start_single` also accepts an explicit component key; its original overload defaults
to the scheduler profile. Library callers that own MPI should initialize it before
starting work. The printing helper does not change the calculation's execution mode.

CTest compares each shipped simulation/diagnostic entry point's information output
with the catalog, checks help/license separation and absence of output files, and
repeats the contract with two ranks when an MPI launcher is available. Startup tests
cover a two-task looper calculation, the embedded scheduler, parallel stdin workers,
and both caller-owned MPI and unconditional cleanup after an information query.

## CFF schema attribution

`cff-1.2.0.schema.json` is an unmodified copy from the Citation File Format project:
https://github.com/citation-file-format/citation-file-format/blob/1.2.0/schema.json

It is distributed under Creative Commons Attribution 4.0 International; see
`CFF-LICENSE`. ALPS's policy schema and generator are separate project code.
