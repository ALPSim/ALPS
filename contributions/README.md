# Contributing an algorithm

Open a pull request that adds your code to `contributions/submission/`: one or
more plain `.cpp` files and a `manifest.yml`. The code needs no knowledge of
this harness; the manifest says which function to call and what it computes,
and the harness writes the glue.

The **Contribution judge** check builds your code, runs it on the test
problems your manifest declares, and passes if at least 50% of them are
answered correctly. Once the PR is merged, **Contribution filing** re-judges it
and moves it to `contributions/verified/<name>_<your GitHub handle>/`, together
with the generated `algorithm.hpp` and a `results.md` with pass rates and
timings by geometry, by quantity and by case. Submitting again under the same
name replaces your own folder.

A PR that adds a submission may not change anything else.

## manifest.yml

```yaml
name: jacobi_ed                  # letters, digits, '_' or '-'
description: One line about the method.
contract: dense_eigenvalues      # how the harness calls your function, see below
entry: jacobi_eigenvalues        # the function to call
models: [tight-binding]          # required
geometry: [1d, 2d]               # optional, default: all
lattices: [open chain lattice]   # optional, default: all
quantities: [energy]             # optional, default: everything the contract provides
inputs: [L, t, N]                # contract model_parameters only
```

Cases outside what the manifest declares are reported as N/A and do not count
for or against the pass rate. Declaring something your code cannot do simply
fails those cases.

| Field | Allowed values |
|---|---|
| `models` | `tight-binding` |
| `geometry` | `1d`, `2d` |
| `lattices` | `open chain lattice`, `chain lattice`, `dimer`, `open square lattice`, `square lattice` (ALPS names) |
| `quantities` | `energy`, `correlation` (G(i,j) = <c+(i) c(j)> in the ground state) |
| `inputs` | `L` (int), `N` (int, particles), `t` (double, hopping), `V` (double, staggered on-site +-V) |

## Contracts

Your function must have exactly this signature.

| `contract` | Signature | Quantities |
|---|---|---|
| `dense_eigenvalues` | `std::vector<double> f(std::vector<double> h, int n)`: all eigenvalues of the symmetric n x n row-major matrix `h`, any order | energy |
| `dense_eigensystem` | `void f(std::vector<double> h, int n, std::vector<double>& values, std::vector<double>& vectors)`: eigenvalues, and normalized eigenvector k at `vectors[k*n .. k*n+n)` | energy, correlation |
| `model_parameters` | `double f(...)`: the ground-state energy, taking the `inputs` in the order listed, with the types above | energy |
| `custom` | none: supply your own `algorithm.hpp` that defines `onboard::makeAlgorithm` returning a `Solver` (see `harness/solver.hpp`) | any |

`contributions/examples/` has one working submission per contract.

## Trying it locally

```bash
bash contributions/harness/judge.sh /tmp/judge path/to/your/submission
```
