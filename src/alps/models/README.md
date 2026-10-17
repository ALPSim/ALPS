# `ALPS::models` — programmatic lattices and models

`<alps/models.hpp>` builds a many-body model in code and compiles it to
numeric local terms. It needs neither the XML catalog nor `alps::Parameters`;
the module depends only on the standard library and Boost headers.
The Python package binds it as `pyalps.models`, with the lattices under
`pyalps.models.lattice` ([pyalps README](../../../python/pyalps/README.md#models-exact-diagonalization-and-classic-dmrg)).

```cpp
#include <alps/models.hpp>

namespace lat = alps::lattice;
namespace mdl = alps::models;

auto model = mdl::heisenberg(lat::chain(8), /*J=*/1.0);   // periodic 8-site chain
mdl::CompiledModel compiled = mdl::compile_hamiltonian(model);
// compiled.op: site_dims and terms (coefficient, dense local factors, JW sites)
```

| Header | Content |
| --- | --- |
| `alps/models/lattice.hpp` | `alps::lattice`: `UnitCell`, `Lattice` (tiling, neighbours, momenta, translation generators) and 15 named lattices |
| `alps/models/finite_graph.hpp`, `geometry.hpp` | `FiniteGraph` and `LatticeGeometry`: owned graph and unit-cell descriptions without a tiler |
| `alps/models/site.hpp` | `LocalOp`, `SiteBasis`, `finalize_site`, built-in spin, boson, hard-core boson, spinless and spinful fermion sites |
| `alps/models/terms.hpp` | `Term`, `TermSum`: operator algebra with fermionic swap signs, Hermitian conjugation, canonicalisation and charge analysis; `+=`/`-=` append in place, so accumulating many terms is linear |
| `alps/models/physical_basis.hpp` | `LocalBasis`, `PhysicalBasis`: ordered local spaces, quantum numbers and parities, independent of any Hamiltonian |
| `alps/models/compiled.hpp` | `CompiledOperator` (dense local factors with Jordan–Wigner strings) and `CompiledModel` (operator plus lattice metadata) |
| `alps/models/hamiltonian.hpp` | `Hamiltonian` (lattice + site bases + bond/site/global term templates), `compile_hamiltonian`, `compile_operator`, `validate_hermiticity` |
| `alps/models/builtin.hpp` | 14 named model factories, from `heisenberg` to `bose_hubbard` and `kondo_lattice` |
| `alps/models/lattice_xml.hpp` | `to_latticegraph_xml`, `to_graph_xml`: the ALPS lattice-library XML of a kit lattice (unit cell with vertex coordinates, extent, boundary); `custom_lattice_name`, `catalog_graph_name` |
| `alps/models/qubit_encoding.hpp` | `QubitEncoding`: the qubit register of a `PhysicalBasis` (binary or Gray code words per site), bit/state-index and statevector interchange |
| `alps/models/pauli.hpp` | `PauliOperator`, `PauliExportOptions`, `export_pauli` (sparse Pauli strings of a `CompiledOperator` or `Hamiltonian`), `local_valid_projector`, `invalid_state_penalty` |

Conventions that consumers rely on:

- Site indices follow the lattice: cells enumerate with the last dimension
  fastest and `site = cell * num_vertices + vertex`.
- Dense local matrices are row-major with row = output state, column = input
  state. A term's amplitude is `coefficient * Π_k op_k(out_k, in_k) * Π_{s ∈ jw} F_s(in_s)`,
  with `F` the diagonal parity operator of the site.
- Charge labels are compared within 1e-9; matrix entries are never pruned by
  that tolerance.

## Pauli export

`export_pauli(compile_hamiltonian(model).op, physical_basis(model))` (or
`export_pauli(model)`) expands a compiled operator into Pauli strings on the
register of a `QubitEncoding`:

```cpp
auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 4.0);
mdl::PauliOperator pauli = mdl::export_pauli(model);      // binary encoding, charges checked
// pauli.num_qubits == 12, pauli.terms: "IIIIIIIIIZZI" -> coefficient, ...
auto const encoding = mdl::make_qubit_encoding(model, mdl::QubitEncodingScheme::gray);
mdl::PauliExportOptions observe;                          // an arbitrary observable:
observe.require_hermitian = observe.require_even = observe.require_conservation = false;
auto const c = mdl::export_pauli(model.operator_decomposition(mdl::TermSum::located(model.site_for(2), "c_up", 2)).op,
                                 mdl::physical_basis(model), observe);
```

- Site `s` occupies `site_qubits[s]` consecutive qubits from `site_offsets[s]`;
  site 0 is the least significant. A string lists qubit 0 at its right end.
  `state_index(local_indices)` is the register index of a physical product
  state; `encode_basis_bits`/`decode_basis_bits` serve registers beyond 63 qubits.
- A local index is written as a binary or Gray code word. The unused code
  words of a site whose dimension is not a power of two carry the
  trace-completed local identity, so the physical block is exact and
  `local_valid_projector`/`invalid_state_penalty` separate the rest.
- Fermionic sites need binary dimension-2/4 occupation bases under the strict
  convention; their Jordan–Wigner strings become `Z` strings and
  `fermion_mode_qubits` lists the modes in site-major (down, up) order.
  Spin-1/2 and fermions may mix; higher spins with fermions are rejected.
- Coefficients at or below `tolerance` are dropped and summed into
  `discarded_l1`, an operator-norm bound on the pruning error. Hermiticity,
  even fermionic parity and conservation of the basis charges are required
  unless switched off; `validation_tolerance` separates those checks from
  pruning. `max_terms` and `max_local_matrix_bytes` bound the expansion and
  the per-site dense transform (O(k 4^k) for k qubits, cached per export).

## Legacy XML bridge

`to_latticegraph_xml(lattice, name)` writes a kit lattice as a `<LATTICES>`
document of the ALPS lattice library: the basis, a unit cell whose typed
vertices carry their fractional `<COORDINATE>` and whose edges carry their
cell offsets, and a `<LATTICEGRAPH>` baking the extent and boundary. The
legacy expansion reproduces the kit's sites, types, coordinates and ordered
bonds. `custom_lattice_name` is a content-addressed name for such documents
and `catalog_graph_name` the catalog graph a factory lattice's boundary
selects. The legacy lattice module's `<alps/lattice/kit.hpp>` builds on these
(and links `ALPS::models`): `resolve_lattice` selects the catalog graph only
when the catalog has it and expands it to the identical structure (the typed
honeycomb and kagome cells, the rectangular and non-cubic 3-D lattices and
custom cells take the custom document), `make_graph_helper` builds the legacy
`graph_helper` in process, `write_lattice_library` writes the document for an
XML application's `LATTICE_LIBRARY`, and `CatalogModel(lattice, parameters)`
runs an XML catalog model on a kit lattice.

## Origin

The code is a port of the lattice and model kit of the experimental
[skilledwolf/alpsx](https://github.com/skilledwolf/alpsx) fork (private at the
time of porting); each file records its source path and the adaptations made.
The fixtures under `tests/fixtures/` were generated by that fork's Python
implementation of the same kit and are frozen here as reference data; the
tests compare the C++ kit against them.

Not ported: the TOML catalog readers and emitters, the bridge to legacy
`alps::Parameters`, the in-process lattice registry and the MPO compiler. The
Pauli/qubit export (`qubit_encoding.hpp`, `pauli.hpp`) follows the fork's,
taking compiled operators instead of raw terms; the XML emitters
(`lattice_xml.hpp`) follow the fork's last XML revision and add the vertex
coordinates, with the parameter bridge in the legacy lattice module.
