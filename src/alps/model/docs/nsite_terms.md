# n-site terms

Terms acting on more than two sites. Supported by exact diagonalization
(`fulldiag`, `sparsediag`). Other applications refuse a model that has
n-site terms on a lattice that has n-site groups.

## Lattice

An `<NSITE>` is an ordered group of distinct sites. It is stored next to the
graph; it does not add bonds.

In a unit cell, each `<SITE>` is a vertex plus a cell offset, as for the
`SOURCE`/`TARGET` of an `<EDGE>`. N-site groups with a site off an open boundary
are dropped.

```xml
<UNITCELL name="ladder cell" dimension="1" vertices="2">
  ...
  <NSITE type="0">
    <SITE vertex="1" offset="0"/>   <!-- role i -->
    <SITE vertex="2" offset="0"/>   <!-- role j -->
    <SITE vertex="2" offset="1"/>   <!-- role k -->
    <SITE vertex="1" offset="1"/>   <!-- role l -->
  </NSITE>
</UNITCELL>
```

In a `<GRAPH>`, sites are vertex ids (1-based):

```xml
<NSITE type="3"><SITE id="1"/><SITE id="4"/><SITE id="2"/></NSITE>
```

## Model

The `sites` attribute names the roles. They are bound to the n-site group's sites
in order.

```xml
<NSITETERM type="0" sites="i j k l">
  K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))
</NSITETERM>
```

- The term is any operator expression over the site basis operators and the
  bond/site operators of the library, like a `BONDTERM`. The number of roles
  must match the n-site group size. N-site groups may have 2, 3, 4 or more sites.
- Site operators may not be raised to a power or used as the argument of a
  function: write `n(i)*n(i)` rather than `n(i)^2`. Such terms are refused
  with an error.
- Without a `type` attribute the term applies to every n-site group type, and `#`
  is replaced by the type, as for site and bond terms. An n-site group type with no
  matching term contributes nothing.
- Fermion signs follow the operator order in the term string and the ALPS
  site order (Jordan-Wigner over sites in index order, up before down on a
  site), the same convention as bond terms.
- Terms need not conserve particle number. If a term maps a basis state out of
  the sector set by `CONSERVED_QUANTUMNUMBERS`, the run stops with an error
  instead of dropping the matrix element.

## Bosons and cutoffs

Boson site bases work as well; there is no sign. The local range is the site
basis's own `0 <= N <= Nmax`. `Nmax` must be set: with its default of
`infinity` the basis cannot be built and the run stops.

Each operator is truncated: `bdag` on `N = Nmax` gives zero, so a matrix
element that would take a site above `Nmax`, even in an intermediate step, is
dropped. For example `b(i)*bdag(i)` on `N = Nmax` gives 0, not `Nmax+1`.

If an n-site term changes a quantum number whose range is set by such a
cutoff parameter (e.g. `K*(bdag(i)*bdag(k) + b(k)*b(i))` changes `N`), the
results depend on the cutoff, and ALPS prints a warning once:

```
Warning: an <NSITETERM> changes the total N, whose local range is cut off at
N <= Nmax = 2. ... check convergence in Nmax.
```

Leave `N` out of `CONSERVED_QUANTUMNUMBERS` for such terms. Fixed ranges
(fermions, hard-core bosons) and ranges set by another quantum number (`Sz` in
`-S..S`) are not cutoffs and give no warning.

## Export to FCIDUMP

`alps2fcidump` writes a fermion model (one spin-1/2 orbital per site, e.g.
`BASIS "fermion"`) with its site, bond and n-site terms as an FCIDUMP file,
the integral format read by quantum chemistry ED, CI and DMRG codes (PySCF,
Molpro, Block2, Dice, NECI, ...):

```
alps2fcidump parm FCIDUMP
```

`parm` is an ordinary ALPS parameter file (`LATTICE`/`GRAPH`, `MODEL`,
couplings, libraries). `Nup_total` and `Ndown_total` give `NELEC` and `MS2`,
or `N_total` and `Sz_total` with a basis such as `"alternative fermion"`;
without them, half filling is written. With `TOTAL_MOMENTUM` the file is still
the real-space Hamiltonian. The export sets `Nup_total`, `Ndown_total`,
`N_total` and `Sz_total` for each sector it builds, so couplings must not
refer to them: `mu="(Nup_total+Ndown_total)/4"` would be exported as
-N^2/4, which equals the model only at that filling.

The same export is a parameter of `fulldiag` and `sparsediag`:
`FCIDUMP="file"` writes the model to `file` before it is diagonalized, with
the same checks. A model that FCIDUMP cannot hold stops the run with the
reason. Each task of a job needs its own `FCIDUMP`; a job in which two tasks
name the same file fails before any task runs. A task that sets only one of
`Nup_total` and `Ndown_total`, to scan the other's sectors, gets a note and
the `NELEC` and `MS2` of the other parameters or of half filling; `alps2fcidump`
refuses it.

The integrals are read off the matrix ALPS builds, in the sectors with at most
one up and one down electron, so they follow ALPS's operator definitions and
fermion signs. Orbital `p` is site `p` (ALPS order, 1-based). A constant term
of the model becomes `E_core`. ORBSYM is all 1.

FCIDUMP holds only a spin-independent two-body Hamiltonian whose (pq|rs) has
the 8-fold symmetry of real orbitals. The tool rebuilds every sector with up
to four electrons (`--check N`) from the integrals and compares it with
ALPS's matrix; if they differ, or a symmetry fails, it writes nothing and
says why. Refused, for example: spin models, a Zeeman field, terms such as the
walkthrough's `K` term that act on up and down electrons differently,
three-body terms, and singlet pair hopping between bonds, which is
spin-independent but has (pq|rs) != (qp|rs).

The check builds ALPS's matrix in the sectors with four electrons, about
n^4/4 states for n orbitals, which dominates the export beyond about 30
orbitals. `--check 3` still catches three-body terms; `--check 2` checks only
the sectors the integrals are read from.

A worked example (one 2x2 plaquette followed through the local table and the
matrix fill), the benchmarks against QDK CASCI and exact results, and the
scripts that produced them are in [`n_body_terms/`](n_body_terms/README.md).
