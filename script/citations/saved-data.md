# Citations in saved ALPS data

ALPS HDF5 checkpoints and results carry the citation recommendations from the
build that produced them. `CITATION.cff` and `CITATIONS.yaml` remain the authority
for new builds. Saved metadata is a historical snapshot, so installing a newer
ALPS version does not change what an old result recommends.

## Stored representation

Metadata lives at `/provenance/alps/citations`, separate from parameters,
observables, random-number state, and restart state. Schema version 1 contains:

```text
/provenance/alps/citations/
  schema_version = 1
  records/<snapshot_sha256>/
    component
    software_version
    catalog_sha256
    selection_basis = "application"
    activity = "calculation" | "analysis" | "unspecified"
    bibliography_cff
    algorithm[]
    implementation[]
    framework[]
    request
    note
    license_note
    notice
    complete = 1
```

`bibliography_cff` is a self-contained, schema-valid CFF 1.2.0 document. It
contains the software metadata, preferred framework citation, and only the
additional papers selected for that component, including declared dependencies.
The role arrays refer to the ordinary ALPS reference identifiers in that CFF.
No custom policy fields are inserted into CFF. This uses CFF's standard
[`preferred-citation` and `references`](https://github.com/citation-file-format/citation-file-format/blob/1.2.0/schema-guide.md).

`notice` is exactly the text produced by that build's `--citations` query. The
separate request and role arrays allow tools to present the same policy without
parsing the notice. All text datasets declare UTF-8 encoding; role keys use
ASCII strings. Empty arrays use the existing ALPS HDF5 null-dataspace convention.

`catalog_sha256` identifies the complete source bibliography and policy, using
SHA-256 of their canonical JSON representation. A record's group name is the
SHA-256 of the record's fields, excluding the group name and completion marker.
Canonical JSON sorts object keys, preserves array order, uses compact separators,
and encodes Unicode as UTF-8 without ASCII escapes. The CFF and notice strings
are included verbatim. This distinguishes changes in policy, bibliography,
software version, component, or activity. It does not identify a source commit
or validate the scientific data; full execution/build provenance is a separate
concern.

The collection is a set of recommendations, sorted by identity on read. It is
not an ordered workflow log. `selection_basis` states what ALPS knows: the
application that produced the contribution. It does not infer a particular
algorithm from simulation parameters. Generic library callers default to the
framework profile with `activity = "unspecified"`; applications explicitly
provide their component and activity.

## Saving, restarting, and deriving results

Repeated saves deduplicate identical records. Loading a checkpoint retains its
records in memory. Saving the resumed result includes those records and the
current build's recommendations. If the catalog or software version changed,
both snapshots remain available; the old one is not rewritten.

Ancestry follows loaded or explicitly supplied data. A fresh result that reuses
an output filename does not inherit that filename's old recommendations. This
requires explicit metadata replacement because ALPS's HDF5 `"w"` mode starts
from a copy of an existing archive. `replace_citations` replaces the set supplied
by a complete-result writer; `write_citations` appends recommendations to an
existing contribution or shared archive. Both reject unsupported schemas before
changing metadata.

The conventional scheduler propagates its application component to task and
worker checkpoints. Parapack clone checkpoints preserve loaded histories and
its evaluator collects clone histories into aggregate HDF5 results. NGS
`mcbase` retains loaded/inherited histories, and application result writers pass
that history to `save_results`. CT-INT, CT-HYB, and standalone Hirsch-Fye carry
input histories into their outputs. DMFT collects its native solver's profile
or citation snapshots supplied in an external solver's returned HDF5 file
before that temporary file is removed. An external solver that supplies no
metadata contributes no inferred profile. MaxEnt retains its HDF5 input's
history and adds a framework analysis snapshot.

Citation writes run inside the existing file-owner/lock boundary. They do not
initialize MPI or change rank ownership. Each separately written worker/clone
checkpoint receives metadata; the existing master writer owns shared results.
Records are marked complete after their fields are written. Readers skip an
unfinished append, and a subsequent append of that identity replaces the partial
record. This marker is not a filesystem transaction or a guarantee against
arbitrary HDF5 file corruption.

Files without metadata remain readable and return an empty citation collection.
Readers do not invent historical recommendations for them. Unsupported schema
versions and structurally malformed metadata raise an error. The schema version
must be a scalar integer, and completion markers must be integer zero or one. The Python reader
also recomputes each snapshot fingerprint; the native reader checks structure,
identifier format, and conflicts when merging records.

## Reading and exporting

```python
import pyalps

records = pyalps.read_citations("result.h5")
for record in records:
    print(record["component"], record["software_version"], record["activity"])
    print(record["notice"])

pyalps.export_citations("result.h5", "CITATION.cff")
```

The reader also accepts an open `pyalps.hdf5.archive` or `h5py.File`. Filename
access uses the existing ALPS bindings; h5py is optional and is not a new pyalps
runtime dependency. Reading needs no YAML parser, catalog lookup, network access,
or write access to the source file.

Export writes the saved CFF document exactly and refuses to overwrite an
existing destination. If a file contains different bibliographies, select a
record explicitly:

```python
pyalps.export_citations("result.h5", "historical.cff", snapshot_id=records[0]["id"])
```

Identical bibliographies shared by calculation and analysis snapshots can be
exported without choosing an identity. Different versions or components are not
silently combined into a new bibliography. To obtain current recommendations,
query the current application's `--citations` explicitly; reading historical
data never refreshes its policy.

Native integrations use `alps/utility/citation_provenance.hpp`:

```cpp
simulation.set_citation_component("interaction");
simulation.inherit_citations(alps::read_citations(input_file));
alps::save_results(results, parameters, output_file,
                   "/simulation/results", simulation.citations());
```

Low-level archive callers use `read_citations`, `merge_citations`,
`make_citation_snapshot`, and `write_citations`/`replace_citations` within their
normal HDF5 save boundary. Pure parameter conversion does not claim a calculation
or add an algorithm citation. Raw archive operations and arbitrary third-party
save implementations must opt in; the archive itself does not guess a component.

## Compatibility and verification

This implementation covers the HDF5 save paths above. XML, XDR, text plots, and
DMRG's private binary scratch files are unchanged; no duplicate sidecar policy
is introduced. Adding metadata to those formats needs a separate compatibility
design. Data copied through a tool that discards `/provenance` also loses its
history, so downstream writers must preserve or explicitly inherit it.

Tests cover generation/CFF validity, old and new build identities, historical
reader behaviour, export selection, malformed and future schemas, incomplete
appends, deduplication, task/worker and NGS checkpoint replacement, fresh output
reuse, clone relocation with serial and two-rank MPI writers, aggregate ancestry
and unchanged observable values, and both external
DMFT solver protocols. Where h5py is available, an
independent reader verifies native-written HDF5 payloads and UTF-8 types against
the generator. Existing serial and MPI CLI tests continue to check the printing
policy independently of storage.
