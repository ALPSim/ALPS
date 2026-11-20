# DCA

`dca` runs the self-consistency loop of the dynamical cluster approximation
(DCA) for Hubbard-type models on clusters from `libcluster/cluster.xml`
(two-dimensional Betts clusters), `cluster1d.xml` and `cluster3d.xml`. It
supports the normal state, the antiferromagnetic and charge-ordered double unit
cell frameworks, and superconducting order in Nambu space. A one-site cluster
gives single-site DMFT on the lattice dispersion.

The impurity problem is solved by an external solver program. `dca` does not
include a solver.

## Provenance

This code is the DCA framework of [CQMP/DMFT](https://github.com/CQMP/DMFT),
imported with its history from branch `master` at commit `642afdf` (June 2026)
by `git subtree`. The upstream code built against ALPSCore. Within ALPS it
builds against `ALPS::params`, `ALPS::gf` and `ALPS::hdf5`, which are copied
from ALPSCore. The following changes were made on import:

- The single-site DMFT loop and its Hilbert transforms (Bethe lattice, density
  of states, two-dimensional band structure), the CDMFT loop and the `SC_LOOP`
  parameter that chose between them were removed.
- The DCA+ loop, auxiliary programs and a test that were not built were removed.
- The code was relicensed under MIT.

## Building

`dca` is built with the ALPS applications when CMake finds Eigen, which
`ALPS::gf` needs. The self-energy interpolation tools `sigma_interpolate_normal`
and `sigma_interpolate_sc` are also built when CMake finds GSL. The cluster
definitions are installed in `share/alps/dca`.

## Running

```sh
dca parameters.ini
```

`dca --help` lists all parameters, `dca --citations` the papers to cite and
`dca --license` the license. The main parameters are:

| Parameter | Meaning |
|---|---|
| `BETA`, `MU`, `U`, `t`, `tprime` | Inverse temperature, chemical potential, interaction, hoppings |
| `N`, `NMATSUBARA`, `FLAVORS` | Imaginary-time points, Matsubara frequencies, flavors (default 2) |
| `dca.LATTICE`, `dca.LATTICE_LIBRARY`, `dca.SITES` | Cluster name, cluster definition file, number of cluster sites |
| `dca.L` | Linear number of lattice momenta per cluster momentum in the coarse-graining |
| `dca.AFM_SELFCONSISTENCY`, `dca.CO_SELFCONSISTENCY`, `sc.NAMBU` | Double unit cell (AFM or charge order) or Nambu frameworks |
| `dca.G0_OMEGA_INPUT`, `dca.G0_OMEGA_OUTPUT` | Starting G0, and G0 written after each iteration to restart from |
| `sc.MAX_IT`, `sc.SELFCONSISTENCY_CONVERGENCE` | Number of iterations, and the change of G at which to stop earlier |
| `sc.WRITE_ITERATIONS` | Write the Green's functions of every iteration (default `true`) |
| `sc.ADJUST_MU`, `sc.TARGET_DENSITY` | Adjust MU to reach a target density |
| `solver.SOLVER` | Command that starts the impurity solver |

## Impurity solver interface

In each iteration, `dca` writes the parameters at `/parameters`, the cluster
G0(k, iωn) at `/G0` and G0(k, τ) at `/G0_tau` to the HDF5 file
`solver.INFILE_H5GF` (default `alps_solver_in.h5gf`). It then runs

```sh
<solver.SOLVER> <parameter file> --MU=<mu>
```

and reads G(k, iωn) at `/G_omega` and G(k, τ) at `/G_tau` from
`solver.OUTFILE_H5GF` (default `alps_solver_out.h5gf`). Both files are removed
after each iteration. The solver reads the same parameter file as `dca`.

Each iteration also writes the Green's functions and self-energy as text
(`G0_omega_<n>`, `G_omega_<n>`, `selfenergy_<n>`, ...) and to `Greens_<n>.h5`,
unless `sc.WRITE_ITERATIONS` is `false`. dca prints the change of G in every
iteration; it stops when the change is below `sc.SELFCONSISTENCY_CONVERGENCE`
or after `sc.MAX_IT` iterations, which is the usual way to run it: run a fixed
number of iterations and judge convergence afterwards. dca says which of the
two ended the run. To continue a run, start dca again with `dca.G0_OMEGA_INPUT` set to the
previous `dca.G0_OMEGA_OUTPUT`; dca refuses a G0 file computed for another
`BETA`, `NMATSUBARA`, cluster or `FLAVORS`.

## Citations

`dca` asks users to cite the review of quantum cluster theories by Maier et
al., Rev. Mod. Phys. 77, 1027 (2005), for the method. There is no paper on this implementation;
its normal-state framework was used in Gull et al., Phys. Rev. B 80, 245102
(2009) and its superconducting framework in Gull, Parcollet and Millis,
Phys. Rev. Lett. 110, 216405 (2013), which are cited for it. Users should
also cite their impurity solver.

## Tests

The cluster geometry and transform tests and an end-to-end test of the loop
run in the `dmft` test component. The end-to-end test runs `dca` at U=0 with a
solver that returns G = G0 on one-site 2D and 3D clusters and on the 4-, 8- and
16-site Betts clusters. It checks the isolated-cluster starting point, a zero
self-energy, convergence at the coarse-grained lattice Green's function, and
the local Green's function against an independent momentum integral. The same
checks run for the double unit cell (AFM and CO, on the bipartite 4-, 8- and
16-site clusters) and Nambu (d-wave pairing field, 1 to 16 sites) frameworks,
with the 2x2 lattice problem integrated independently. A second end-to-end
test runs `dca` with a solver that returns G = [G0^-1 - Sigma]^-1 for a fixed
Sigma(iw) = S0 + S1/(iw) in all four frameworks (with d-wave pairing in Nambu
space). dca's self-energy must be that Sigma in every iteration, and the local
Green's function that of the lattice with Sigma. The solver also checks that
the G0(tau) it receives is the transform of G0(iw), and with a Sigma that
differs between the spins, that dca hands the spin average of G to Dyson's
equation. The bare
tail test checks the high-frequency tails attached to G0 against G0 itself in
the normal-state, double unit cell and Nambu frameworks.

```sh
ctest -L dmft -R dca_
```
