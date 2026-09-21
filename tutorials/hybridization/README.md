# Hybridization-expansion solver reference

The CT-HYB solver computes quantum impurity models using continuous-time hybridization-expansion Monte Carlo. Start with the [Python interface tutorial](../hybridization-01-python/), then continue with [Kondo physics](../hybridization-02-kondo/), [retarded interactions](../hybridization-03-retarded-interaction/), and [multiorbital spin freezing](../hybridization-04-spinfreezing/).

## Build and run

Follow the [ALPS build instructions](../../CONTRIBUTING.md#getting-started-with-the-code) with `ALPS_BUILD_APPLICATIONS=ON`. The native executable is named `hybridization`; `p2h5` converts its text parameters to HDF5. The [pyalps build instructions](../../python/pyalps/README.md) explain how to build the Python interface against an installed SDK.

For a plain build directory:

```sh
cmake --build build --target hybridization p2h5
```

After installation, with the SDK's `bin` directory on `PATH`, run a prepared input as follows:

```sh
p2h5 hyb1.h5 < hyb1.param
hybridization hyb1.h5
```

Run Python examples from their tutorial directory using an environment containing pyalps. Serial calls to `pyalps.cthyb.solve(parameters)` require no MPI import. MPI support for native runs is enabled explicitly when building ALPS.

## Detailed manual

The [LaTeX manual](hybdoc.tex) documents hybridization and retarded-interaction inputs, parameters, measurements, self-consistency, analysis, and the solver's limitations. Its [bibliography](refs.bib) accompanies it. This scientific reference is maintained beside the tutorials; generated PDFs are not stored in the source tree.

To render the manual, install a TeX distribution with REVTeX 4.1 and run these commands in this directory:

```sh
pdflatex hybdoc.tex
bibtex hybdoc
pdflatex hybdoc.tex
pdflatex hybdoc.tex
```
