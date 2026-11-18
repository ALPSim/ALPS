# Result schema for Monte Carlo observables

Status: `alps::mcbase` simulations, the CT-HYB and CT-INT solvers and
`pyalps.ngs` write this schema. The scheduler and parapack applications still
write the layout of 3.0 and follow. The conformance tests in `tests/statistics`
are written against it.

ALPS has three implementations of Monte Carlo statistics today: the legacy
classes behind `<alps/alea.h>`, `alps::accumulators`, and this directory. They
store results in three different HDF5 layouts. The plan is to keep one
implementation (this one, with an ALPS observable layer on top) and one layout.
This file defines that layout and the rules the implementation follows.

## Design rules

1. **Safe defaults.** The default error of an autocorrelated observable is
   taken from the longest bin length that still has at least 128 bins. This is
   the rule the legacy classes and `alps::accumulators` use. If no bin length
   above 1 qualifies, the autocorrelation time is NaN and the result is marked
   as not converged. A naive error is never reported with tau = 0.
2. **Equal-length, time-ordered bins.** An observable keeps up to B bins of
   equal length plus one partial bin. When the bins are full, neighbours are
   merged and the bin length doubles. Merging runs concatenates their bins.
3. **One default per observable:** count, mean, error, tau and convergence flag
   from a binning analysis, plus 128 bins for jackknife. An observable can be
   declared without bins.
4. **Textbook estimators.** The standard error from M bin means is
   sqrt(sum of squared deviations / (M (M - 1))). The legacy classes multiply
   this by sqrt(N / (N - 1)) for N samples; that factor is not kept.
5. **Functions of results are jackknife over aligned bins, or an error.**
   `a * b`, `a / b` and `f(a)` need bins on every operand, with equal counts and
   lengths (after re-binning to the coarser length). There is no fallback that
   treats operands as uncorrelated.
6. **Results and checkpoints are different things on disk.** Results follow this
   schema. Checkpoints hold the state needed to continue a run and carry no
   compatibility promise between minor releases.
7. **One vocabulary.** The legacy class names are not kept as aliases.
8. **Sizes are learned at the first sample,** as with legacy vector
   observables, unless the declaration fixes them.

## Results

One group per observable, in the places results are stored today, for example
`/simulation/results/<name>` and
`/simulation/realizations/0/clones/<n>/results/<name>`.

| Path | Compared with 3.0 | Content |
| --- | --- | --- |
| `@alea_schema` | new | integer, `1` |
| `@kind` | new | `scalar`, `vector` or `histogram` |
| `@value_type` | new | `real` or `complex` |
| `count` | unchanged | number of samples |
| `mean/value`, `mean/error` | unchanged | scalar, or 1-d array for a vector |
| `mean/error_convergence` | unchanged | 0 converged, 1 maybe converged, 2 not converged |
| `variance/value` | unchanged | variance of a single sample |
| `tau/value` | unchanged | integrated autocorrelation time; NaN if not available |
| `labels` | unchanged | labels of vector components |
| `timeseries/data`, `@binsize`, `@binningtype` | see below | equal-length, time-ordered bin means |
| `@sign`, and the sibling `"<sign> * <name>"` | unchanged | signed observables |
| `histogram`, `@min`, `@max`, `@stepsize` | unchanged | integer histograms |
| `binning/binsize`, `binning/error` | new | error of the mean at each bin length |
| `covariance/value` | new, optional | vector observables that ask for it |

`timeseries/data` holds bin means. In 3.0 files it has two meanings: evaluated
results (also those of ALPSCore) hold bin means there and either mark the
group with `@cannotrebin` or do not record the bin length, while the
observables a simulation writes while running hold bin sums together with
`@binsize`. The reader tells the cases apart by `@alea_schema`,
`@cannotrebin` and `@binsize`, as ALEA's `mcdata` did.

No longer written, because they are internals of one implementation: `sum`,
`sum2`, `timeseries/logbinning`, `timeseries/logbinning2`,
`timeseries/logbinning_lastbin`, `timeseries/logbinning_counts`,
`timeseries/partialbin`, `timeseries/partialbin2`, `timeseries/data2`.
`binning/*` carries what `pyalps.loadBinningAnalysis` computes from them.

The convergence flag follows the legacy classes (`SimpleBinning::converged_errors`).
Let e be the error at the chosen bin length. The errors at the three next
shorter bin lengths are looked at in turn, shortest first: one that is at least
e sets the flag to converged, one below 0.824 e sets it to not converged, and
one below 0.9 e sets it to maybe converged unless it is already not converged.
With fewer than four bin lengths the flag is maybe converged.

## Checkpoints

Restart state is stored beside the results, not inside them:

- `/simulation/realizations/0/clones/<n>/checkpoint/<name>` for scheduler and
  parapack clones,
- `/simulation/realizations/0/clones/0/measurements/<name>` for `alps::mcbase`
  simulations.

Each group carries `@alea_state` (an integer version) and the accumulator type.
A checkpoint is read only by the minor version that wrote it; anything else is
refused with a message that says so.

## Readers

- Schema 1 files.
- Files written by ALPS 2.x and 3.0: the kind is inferred from the datasets that
  are present, as the legacy reader does, and `mean/value`, `mean/error`,
  `tau/value` and the convergence flag are returned as stored.
  `tests/statistics/fixtures/legacy_observables.h5` is such a file, with the
  numbers a reader must return in `legacy_observables.expected`.
- `pyalps.loadMeasurements` needs `count`, `mean/*` and `timeseries/data`, which
  do not change. `pyalps.loadBinningAnalysis` reads `binning/*` from schema 1
  files and keeps its existing code for older files.

## Changes a user can see

- Error bars change by the relative amount 1/(2N) for N samples (rule 4).
- A short run reports NaN for tau and "not converged" where the unrepaired code
  in this directory reports the naive error (rule 1).
- Scripts that read `sum`, `sum2` or `timeseries/logbinning*` directly from a
  result file need to use `binning/*` or the pyalps loaders.
