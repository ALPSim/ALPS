# Command-line tools

These directories own executable consumers of ALPS modules. They do not add
new independently linkable libraries. The root CMake file dispatches to each
group; executable names and installation components remain stable.

| Group | Built or installed programs | Preserved inactive sources |
| --- | --- | --- |
| `parameters/` | `parameter2xml`, `parameter2hdf5`, `p2h5` | None |
| `lattice/` | `lattice2xml`, `printgraph` | `pltgraph.py` and historical local fixtures |
| `scheduler/` | `convert2xml`, `compactrun`, `snap2vtk` | None |
| `parapack/` | `pevaluate`, `poutput` | `xml2archive.C` |
| `diagnostics/` | `pconfig` | None |
| `xml/` | `alps-xml` on Unix | `txt2archive.C` and all historical shell wrappers |
| `result_archive/` | None | SQLite-backed XML result indexing and plotting sources |
| `alea/` | None | C++ and Python mean/variance analysis programs |
| `launchers/` | None | Historical `alpspython` shell and Windows templates |
| `installer/` | None | Historical macOS postflight template |

All eleven C++ commands install in the `tools` component; `alps-xml` installs in
the `xml` component. `pconfig` links `ALPS::utilities`; the other C++ commands
continue to use `ALPS::alps`. Source ownership is registered for all C++ groups,
including inactive files, without enabling additional programs or dependencies.

The parameter tools retain the older `alps::Parameters` grammar, job generation
and seeding behavior. The result archive group is separate from the HDF5
implementation. Inactive sources are preserved for a separate functionality
assessment; their presence does not claim current build or runtime support.
Central CLI integration tests and their fixtures remain under `tests/cli/`.
