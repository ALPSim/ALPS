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
| `xml/` | `alps-xml` on Unix | `txt2archive.C` |
| `installer/` | None | Historical macOS postflight template |

C++ commands install in the `tools` component; `alps-xml` installs in the `xml` component. `pconfig` links `ALPS::utilities`; the other C++ commands use `ALPS::alps`.

The parameter tools retain the older `alps::Parameters` grammar, job generation
and seeding behavior. Remaining inactive sources are preserved for a separate functionality
assessment; their presence does not claim current build or runtime support.
Central CLI integration tests and their fixtures remain under `tests/cli/`.
