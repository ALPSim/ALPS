# Monte Carlo simulation base

This directory holds two implementations of the Monte Carlo simulation base
class `alps::mcbase`, its MPI adapter `alps::mcmpiadapter`, the schedule and
stop callbacks, and `alps::random01`.

- The NGS implementation, ALPS's own: `<alps/mcbase.hpp>`,
  `<alps/mcmpiadapter.hpp>`, `<alps/check_schedule.hpp>`,
  `<alps/stop_callback.hpp>`, `<alps/progress_callback.hpp>` and
  `<alps/ngs/api.hpp>`, with sources in `src/ngs/`. Its observables are
  ALEA's `alps::mcobservables`. It is compiled into `ALPS::alps`.
- `ALPS::mc`, the Monte Carlo module of the
  [ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore):
  `<alps/mc/*.hpp>`, with sources in `src/` and tests in `tests/`. Its
  observables are ALPSCore's accumulators. It is copied from ALPSCore at commit
  `7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25).

Both define the same classes in namespace `alps`, so a program uses one or the
other. For now `ALPS::mc` is built and tested in the source tree only; it is
not part of `ALPS::alps` or the installed SDK.

## Provenance

The imported revision is under the MIT license in [LICENSE.TXT](LICENSE.TXT).
Every imported file refers to [COPYRIGHT.TXT](COPYRIGHT.TXT),
[LICENSE.TXT](LICENSE.TXT) and [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT), which are
ALPSCore's top-level files, copied unchanged. Keep the notices when modifying
these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `mc/include/alps/mc/*.hpp` (6 files) | `include/alps/mc/` | Copied unchanged |
| `mc/src/*.cpp` (3 files) | `src/` | Copied unchanged |
| `mc/test/*.cpp` (9 files) | `tests/` | Copied unchanged |
| `COPYRIGHT.TXT`, `LICENSE.TXT`, `ACKNOWLEDGE.TXT` | here | Copied unchanged |

Not copied: ALPSCore's build files and `alps-mc.pc.in`. The CMake files and
this README are ALPS's own.

The tests run as one CTest test per GoogleTest case, and the MPI tests on 2
and 3 ranks. `reduce` and `signed_obs` have no GoogleTest cases; their exit
status is the check.
