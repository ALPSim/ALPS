# Monte Carlo simulation base

`ALPS::mc` is the Monte Carlo module of the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore): the simulation
base class `alps::mcbase`, its MPI adapter `alps::mcmpiadapter`,
`alps::check_schedule`, `alps::stop_callback`, `alps::random01` and the
functions of `<alps/mc/api.hpp>`, in `<alps/mc/*.hpp>`, with sources in `src/`
and tests in `tests/`. Its observables are ALPSCore's accumulators. It is
copied from ALPSCore at commit `7146b9e1f017938a94e5dae35d88467cc5ba7969`
(2026-09-25). `ALPS::alps` links it, and the SDK installs it.

It replaces ALPS's own implementation of the same classes, the NGS one
(`<alps/mcbase.hpp>`, `<alps/mcmpiadapter.hpp>`, `<alps/check_schedule.hpp>`,
`<alps/stop_callback.hpp>`, `<alps/ngs/api.hpp>`, `<alps/ngs/random01.hpp>`),
whose observables were ALEA's `alps::mcobservables`. The
[changelog](../../../CHANGELOG.md) lists what a simulation changes.

ALPS's own files here:

- `<alps/progress_callback.hpp>`, which prints the completed fraction of a
  simulation. It keeps the schedule of the NGS `alps::check_schedule`, which
  could report the time to the next check.
- `tests/runner*.cpp`, tests of a simulation run, checkpointed and resumed, on
  one process and over MPI. They were the NGS scheduler's tests.
- `<alps/ngs/detail/export_sim_to_python.hpp>` and
  `<alps/python/save_observable_to_hdf5.hpp>`, which forward to pyalps's
  headers for older Python exports.

## Provenance

The imported revision is under the MIT license in [LICENSE.TXT](LICENSE.TXT).
The imported copyright and license notices are retained in [COPYRIGHT.TXT](COPYRIGHT.TXT) and [LICENSE.TXT](LICENSE.TXT). The SDK installs both under `share/alps/licenses/mc/`. Publication pointers in imported comments have been removed; the source relationships below refer to the remaining code.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `mc/include/alps/mc/*.hpp` (6 files) | `include/alps/mc/` | Copied, except that `mcbase.hpp` gains a virtual destructor (see below) |
| `mc/src/*.cpp` (3 files) | `src/` | Copied |
| `mc/test/*.cpp` (9 files) | `tests/` | Copied |
| `COPYRIGHT.TXT`, `LICENSE.TXT` | here | Copied |

Not copied: ALPSCore's build files and `alps-mc.pc.in`. The CMake files and
this README are ALPS's own.

The tests run as one CTest test per GoogleTest case, and the MPI tests on 2
and 3 ranks. `reduce` and `signed_obs` have no GoogleTest cases; their exit
status is the check. The MPI tests are listed in the MPI job of
`.github/workflows/ci.yml`, which builds only the targets it names.

### Local changes

`alps::mcbase` has virtual member functions but, in ALPSCore, no virtual
destructor. `mcbase.hpp` declares one, marked with a comment, so a simulation
can be destroyed through a pointer to `mcbase`. pyalps relies on this: a
Python simulation is a nanobind trampoline derived from `mcbase`, and nanobind
destroys it as an `mcbase`.
