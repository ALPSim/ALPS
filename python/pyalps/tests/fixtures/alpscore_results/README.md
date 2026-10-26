# ALPSCore results fixture

`results_7146b9e1.h5` was written by
[ALPSCore](https://github.com/ALPSCore/ALPSCore) at commit 7146b9e1, using
`writer.cpp` in this directory. It is laid out as ALPSCore's Monte Carlo
programs write their output: ALPSCore's parameters under `/parameters`, and a
`result_set` of mean, no-binning, log-binning and full-binning accumulators,
scalar and vector, under `/simulation/results`.

The writer then reads the file back with ALPSCore's archive and prints the
stored counts, means, errors and binning analyses; that output is
`expected.json`. `test_alpscore_results.py` checks that
`pyalps.loadMeasurements` and `pyalps.loadBinningAnalysis` return those values.

ALPSCore's results differ from ALEA's where pyalps looks:

| Result | ALEA | ALPSCore |
| --- | --- | --- |
| Binning analysis | `timeseries/logbinning`, `logbinning2` and `logbinning_counts`, from which pyalps computes the errors | `error_bins`: the error at each binning level |
| Autocorrelation time | `tau/value` | `tau` |
| Variance | `variance/value` | not stored |
| Full-binning time series `timeseries/data` | bin sums | bin means, marked with `@cannotrebin` |

Count, `mean/value`, `mean/error` and the time series' `@binsize` and
`@maxbinnum` attributes match. ALEA's `MCScalarData` and `MCVectorData`
divide the bins by `@binsize` unless `@cannotrebin` is present, so they load
ALPSCore's bin means unchanged. They do not see ALPSCore's `tau`.

Keep the file frozen; a copy regenerated with ALPS's own accumulators would no
longer test compatibility with ALPSCore.

## Regenerating

`writer.cpp` is not part of the ALPS build. With ALPSCore's `utilities`,
`hdf5`, `params` and `accumulators` modules installed under `$ALPSCORE`, build
it against ALPSCore and run it here:

```sh
c++ -std=c++17 -I"$ALPSCORE/include" writer.cpp -o writer_core \
  -L"$ALPSCORE/lib" -lalps-accumulators -lalps-params -lalps-hdf5 -lalps-utilities -lhdf5
./writer_core results_7146b9e1.h5 > expected.json
```

The samples come from `std::mt19937`'s raw output, which the C++ standard
fixes, so a rebuild on another platform writes the same values.
