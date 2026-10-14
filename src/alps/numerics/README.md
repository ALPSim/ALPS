# Numerics

`ALPS::numerics` is an interface target with ALPS's numerical helpers, array
mathematics and matrix/vector interfaces under `alps/numeric/`.

## Headers from ALPSCore

Some headers are ALPSCore's, so that ALPSCore code builds against ALPS
unchanged. They are copied from ALPSCore commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969` (`utilities/include/alps/numeric/`).
Their notices refer to ALPSCore's top-level
[`COPYRIGHT.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/COPYRIGHT.TXT),
[`LICENSE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/LICENSE.TXT)
and
[`ACKNOWLEDGE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/ACKNOWLEDGE.TXT).

| Header | Relationship |
| --- | --- |
| `inf.hpp` | Copied unchanged |
| `vector_functions.hpp` | ALPSCore's version, which replaced ALPS's; ALPS additions at the end (see below) |

### `vector_functions.hpp`

ALPS's older version brought Boost's `std::vector` operators into
`alps::numeric` with a using-directive. ALPSCore's declares its own operators
there instead. They treat a default-constructed vector as the zero vector, and
compound assignment throws `std::runtime_error` on a size mismatch. ALPSCore
also moved the `plus`, `minus`, `multiplies`, `divides` and `unary_minus`
functors here from `functional.hpp`, which now only includes this header.

ALPS adds, at the end of the header, what ALPS code used from the older
version: scaling by a scalar of another arithmetic type (`vector<double> /
std::size_t`), which Boost's operators allowed, unary `+`, and `operator<<` for
vectors.

A scope that declares both Boost's and ALPS's vector operators, with
`using boost::numeric::operators::operator*;` and `using
alps::numeric::operator*;`, now finds two equally good candidates for vector
arithmetic and fails to compile. Keep only the `alps::numeric` declarations.
