# Citing ALPS

If ALPS contributes to published research, please cite the framework papers below and any method-specific paper relevant to the application you used. Citation is requested as scholarly acknowledgement; it is not a condition of the MIT license in [`LICENSE.txt`](LICENSE.txt).

## Framework papers

Please cite both ALPS framework papers:

1. A. F. Albuquerque *et al.*, “The ALPS project release 1.3: Open-source software for strongly correlated systems,” *Journal of Magnetism and Magnetic Materials* **310**, 1187–1193 (2007). [doi:10.1016/j.jmmm.2006.10.304](https://doi.org/10.1016/j.jmmm.2006.10.304)
2. B. Bauer *et al.*, “The ALPS project release 2.0: Open source software for strongly correlated systems,” *Journal of Statistical Mechanics: Theory and Experiment* **2011**, P05001 (2011). [doi:10.1088/1742-5468/2011/05/P05001](https://doi.org/10.1088/1742-5468/2011/05/P05001)

## Method-specific papers

Add the applicable paper from this table. Application names match directories or executables in this repository.

| Application or component | Additional citation |
| --- | --- |
| `applications/qmc/looper` (`loop`, `loop_mpi`) | S. Todo and K. Kato, “Cluster Algorithms for General-S Quantum Spin Systems,” *Physical Review Letters* **87**, 047203 (2001). [doi:10.1103/PhysRevLett.87.047203](https://doi.org/10.1103/PhysRevLett.87.047203) |
| `applications/qmc/qwl` | M. Troyer, S. Wessel, and F. Alet, “Flat Histogram Methods for Quantum Systems,” *Physical Review Letters* **90**, 120201 (2003). [doi:10.1103/PhysRevLett.90.120201](https://doi.org/10.1103/PhysRevLett.90.120201) |
| Continuous-time QMC impurity solvers and the DMFT framework | E. Gull, P. Werner, S. Fuchs, B. Surer, T. Pruschke, and M. Troyer, “Continuous-time quantum Monte Carlo impurity solvers,” *Computer Physics Communications* **182**, 1078–1082 (2011). [doi:10.1016/j.cpc.2010.12.050](https://doi.org/10.1016/j.cpc.2010.12.050) |
| DMRG/MPS applications | M. Dolfi *et al.*, “Matrix product state applications for the ALPS project,” *Computer Physics Communications* **185**, 3430–3440 (2014). [doi:10.1016/j.cpc.2014.08.019](https://doi.org/10.1016/j.cpc.2014.08.019) |

The framework papers are sufficient for `sse`, `spinmc`, `fulldiag`, `sparsediag`, and `worm` unless a publication describes a more specific algorithmic reference. For specialized models or algorithms, also cite the primary scientific source on which the calculation is based.

## Citation metadata

Use the DOI links above to obtain current BibTeX, RIS, or other citation metadata from the publishers. This avoids maintaining duplicate hand-written records that can drift from the authoritative metadata.

When describing reproducibility, also record the ALPS release or Git commit used in the calculation.
