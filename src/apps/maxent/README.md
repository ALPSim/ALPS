# MaxEnt

The implementation in `src/` is adapted from
[CQMP/Maxent](https://github.com/CQMP/Maxent) at commit
[`a3992534228b0de732e033ae23154969724814be`](https://github.com/CQMP/Maxent/commit/a3992534228b0de732e033ae23154969724814be).
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
  Application builds additionally require Eigen >= 3.3. The Legendre kernel
  uses the existing Boost.Math Gauss-Kronrod quadrature instead of upstream's
  GSL call, with an equivalent integrand that avoids exponential overflow.
  Library-only builds do not require Eigen.
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
- Size the LAPACK SVD buffers for rectangular kernels and report failed
  decompositions; initialize omitted text covariance entries to zero.
- Keep odd logarithmic grids within the requested frequency endpoints and
  half-Lorentzian grids away from the tangent pole at coarse resolutions.
  Include the first frequency bin in complex kernels, correcting its omission
  in the legacy solver; complex spectra can therefore differ numerically.

The replacement also brings CQMP's fermionic Legendre kernel, additional default
models, and spectrum backcontinuation. Bosonic Legendre input is rejected because
the imported quadrature implements only the fermionic kernel. CQMP's standalone
parameter-file CLI and `MODEL_RUNS` loop are not part of the ALPS entry point.
Scientific output paths and `NORM` scaling are covered by the shared native/Python
references.

## License and scientific credit

CQMP/Maxent's maintainers agreed to MIT licensing for this reintegration, as
confirmed by the reintegration requester on 2026-10-14. `LICENSE.TXT` records
those MIT terms. The pinned upstream tree still contains its earlier GPL text;
the MIT terms here rely on that maintainer agreement. The original source
copyright notices and `ACKNOWLEDGE.TXT` are preserved. This agreement covers
Maxent; external dependencies retain their own licenses.

The citation catalog credits Levy, LeBlanc and Gull's 2017 implementation and
the Bryan and Jarrell–Gubernatis method references. ALPS supplies the current
framework citation; the preserved upstream acknowledgment file also records
its historical ALPS/ALPSCore requests. Run `maxent --citations` for the current
application guidance. Installed SDKs and Python packages include the MaxEnt
license and acknowledgment files.

## Validation

Run `ctest --test-dir <build> -R maxent --output-on-failure` for the public API,
CLI, grid/kernel checks and existing physical reference cases. The Python
binding uses the same reference checks in
`python/pyalps/tests/test_maxent_reference.py`.
