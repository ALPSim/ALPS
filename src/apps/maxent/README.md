# MaxEnt

The implementation in `src/` is adapted from
[CQMP/Maxent](https://github.com/CQMP/Maxent) at commit
[`e727650a1f0a427e95efebc409d907b4c84011ec`](https://github.com/CQMP/Maxent/commit/e727650a1f0a427e95efebc409d907b4c84011ec)
on its `modernize/step2` branch. Upstream `master` is still at the 2018 commit
`a3992534228b0de732e033ae23154969724814be`; the later commits fix kernels,
grids and bootstrap seeding and replace GSL and most Boost utilities.
It replaces the older ALPS MaxEnt implementation. The ALPS executable, Python
`AnalyticContinuation` binding and `alps::solvers::maxent(params, output_file)`
entry point continue to use the same interface.

## Import and local adaptations

All upstream `src/` files except its standalone `maxent.cpp` are imported with
their copyright notices. ALPS retains its own `cli/maxent.cpp` and
`src/maxent_run.cpp`. The companion Padé, Kramers–Kronig and Legendre conversion
programs and the upstream build system are not imported. The grid tests identify
their upstream source separately.

Legacy ALPS portions are retained from commit
`d4533fd30c65656cf5ab7fa448dd6b921b39aa6f`, with their original copyright and MIT
notices alongside the upstream notices. `maxent.hpp`, `maxent_helper.cpp` and
`maxent_simulation.cpp` retain their interfaces or diagnostic/output conventions.
The former `maxent_parms.hpp/.cpp` declarations, input handling, kernels and grids
are now in `maxent_params.hpp/.cpp`, `maxent_kernel.cpp` and `maxent_grid.cpp`;
the former `default_model.hpp` portions are split between `default_model.hpp`
and `default_model.cpp`.

The adaptations are confined to this application:

- Use ALPS's existing typed parameters and HDF5 archive instead of requiring a
  separate ALPSCore installation. Parameter defaults are supplied locally.
- Use Eigen for matrix operations and ALPS's existing LP64 LAPACK for the SVD.
  Application builds additionally require Eigen >= 3.3. Library-only builds do
  not require Eigen.
- Preserve the explicit output filename, stop callbacks and signal setup of
  the ALPS solver. Text output and backcontinuation remain opt-in; set
  `BACKCONTINUE=1` to enable the latter. Backcontinuation with a covariance
  matrix is rejected because the imported routine does not undo the covariance
  basis rotation. The default frequency grid remains linear.
- Preserve ALPS's particle-hole symmetric and complex bosonic frequency
  kernels and `pi * omega * spectrum` bosonic output convention. These differ
  from the pinned CQMP implementation.
- Preserve the indexed three-column text input, uniform imaginary-time grid without
  `TAU_i`, inline-data covariance input, the upper-triangle covariance
  convention, `ENFORCE_NORMALIZATION`, and the `deltaOmega.dat` diagnostic.
- Preserve the default-model column in spectrum text output and linear
  extrapolation beyond tabulated default-model endpoints. For HDF5 input,
  specifying `COVARIANCE_MATRIX`, even with an empty value, selects
  `/Covariance`; omitting it selects `/Error`.
- Evaluate the imaginary-time bosonic kernel at |omega|, where it is even and
  no exponential overflows. Upstream maps the overflow's NaN to zero, which is
  wrong near tau = 0 and beta once beta * |omega| exceeds about 709.
- Seed `GENERATE_ERR` bootstrap errors with `SEED`, default 42 as in ALPS's
  `mcbase`, rather than upstream's 0. The normal variates come from the
  C++ standard library, so seeded errors are reproducible for a given standard
  library, not across platforms.
- Size the LAPACK SVD buffers for rectangular kernels and report failed
  decompositions; initialize omitted text covariance entries to zero.
- Keep odd logarithmic grids within the requested frequency endpoints and
  half-Lorentzian grids away from the tangent pole at coarse resolutions. Both
  defects remain upstream.
  Include the first frequency bin in complex kernels, correcting its omission
  in the legacy solver; complex spectra can therefore differ numerically.

The replacement also brings CQMP's fermionic Legendre kernel, additional default
models, and spectrum backcontinuation. Bosonic Legendre input is rejected because
upstream applies its fermionic Legendre kernel to bosonic data. CQMP's standalone
parameter-file CLI and `MODEL_RUNS` loop are not part of the ALPS entry point.
Scientific output paths and `NORM` scaling are covered by the shared native/Python
references.

## License and scientific credit

The imported code is under the MIT license in [LICENSE.TXT](LICENSE.TXT). Its
copyright notices and `ACKNOWLEDGE.TXT` are preserved.

The citation catalog credits Levy, LeBlanc and Gull's 2017 implementation and
the Bryan and Jarrell–Gubernatis method references. ALPS supplies the current
framework citation; `ACKNOWLEDGE.TXT`, kept from the 2018 import (upstream now
has `CITATION.md` with the same requests), also records its historical
ALPS/ALPSCore requests. Run `maxent --citations` for the current
application guidance. Installed SDKs and Python packages include the MaxEnt
license and acknowledgment files.

## Validation

Run `ctest --test-dir <build> -R maxent --output-on-failure` for the public API,
CLI, grid/kernel checks and existing physical reference cases. The Python
binding uses the same reference checks in
`python/pyalps/tests/test_maxent_reference.py`.
