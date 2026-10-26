# Legacy interfaces

Deprecated ALPS 3.0 C++ spellings, kept so existing downstream code keeps
building while it migrates. Everything here is wired into the build by one
file, `ALPSLegacy.cmake`, behind the `ALPS_LEGACY_API` option (default `ON`).
With the option off, none of these files are on the include path or installed.

Not every interface removed since 3.0.0 has a forwarder here. The NGS
parameters, Monte Carlo base, observables and accumulators are retired
outright, because ALPSCore's replacements use the same class names;
`CHANGELOG.md`'s Migration section records that exception and lists the
retired interfaces.

| Path | What it keeps working | Replacement |
| --- | --- | --- |
| `include/alps/ngs/{cast,stacktrace,signal,stringify,short_print}.hpp`, `include/alps/ngs/detail/{remove_cvr,type_wrapper}.hpp` | the old include spellings | `<alps/utilities/...>` with the same file name |
| `include/alps/utility/{size,temporary_filename,make_copy}.hpp`, `include/alps/utility/{copyright,citations,os,vmusage,vectorio,encode,resize,data}.hpp` | the old include spellings | `<alps/utilities/...>` with the same file name |
| `include/alps/utility/numeric_cast.hpp`, `include/alps/utility/cli.hpp` | the old include spellings | `<alps/numeric/numeric_cast.hpp>` (`ALPS::numerics`), `<alps/cli.hpp>` (`ALPS::cli`) |
| `include/alps/utility/{bitops,generator}.hpp` | the headers themselves, which have no users left | none; `bitops` is superseded by C++20 `<bit>` |
| `include/alps/legacy/ngs_{signal,stacktrace,stringify}.hpp` | `alps::ngs::signal`, `alps::ngs::stacktrace()`, `ALPS_NGS_STRINGIFY` | `alps::signal`, `alps::stacktrace()`, `ALPS_STRINGIFY` |

Each forwarding header prints a `#pragma message` naming its replacement.
`ngs/cast.hpp`, `utility/resize.hpp` and `utility/size.hpp` also include the
Boost headers that their ALPS 3.0 versions included and the canonical headers
dropped, so code relying on those transitive includes keeps building. The
canonical headers include the `alps/legacy/` name headers through a
`__has_include` check, so the old names exist exactly when this layer does;
using them triggers `-Wdeprecated-declarations`. Macros cannot be marked
deprecated, so `ALPS_NGS_STRINGIFY` only disappears with the layer.

## Tests

`tests/` builds the old spellings and names in the build tree (CTest label
`legacy`).

## Removing the legacy layer

Delete this directory. The remaining hooks are inert without it:
`include(legacy/ALPSLegacy.cmake OPTIONAL)` in the top-level `CMakeLists.txt`,
`include(ALPSLegacyConfig.cmake OPTIONAL)` in `cmake/ALPSConfig.cmake.in`, and
the `__has_include(<alps/legacy/...>)` checks in the canonical headers.

This directory follows the layout of the legacy layer proposed in
[upstream pull request 209](https://redirect.github.com/ALPSim/ALPS/pull/209); whichever
change lands second merges into the other's `ALPSLegacy.cmake` and README.
