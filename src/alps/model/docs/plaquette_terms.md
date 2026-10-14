# n-site (plaquette) terms

Terms acting on more than two sites. Supported by exact diagonalization
(`fulldiag`, `sparsediag`). Other applications refuse a model that has
plaquette terms on a lattice that has plaquettes.

## Lattice

A `<PLAQUETTE>` is an ordered group of distinct sites. It is stored next to the
graph; it does not add bonds.

In a unit cell, each `<SITE>` is a vertex plus a cell offset, as for the
`SOURCE`/`TARGET` of an `<EDGE>`. Plaquettes with a site off an open boundary
are dropped.

```xml
<UNITCELL name="ladder cell" dimension="1" vertices="2">
  ...
  <PLAQUETTE type="0">
    <SITE vertex="1" offset="0"/>   <!-- role i -->
    <SITE vertex="2" offset="0"/>   <!-- role j -->
    <SITE vertex="2" offset="1"/>   <!-- role k -->
    <SITE vertex="1" offset="1"/>   <!-- role l -->
  </PLAQUETTE>
</UNITCELL>
```

In a `<GRAPH>`, sites are vertex ids (1-based):

```xml
<PLAQUETTE type="3"><SITE id="1"/><SITE id="4"/><SITE id="2"/></PLAQUETTE>
```

## Model

The `sites` attribute names the roles. They are bound to the plaquette's sites
in order.

```xml
<PLAQUETTETERM type="0" sites="i j k l">
  K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))
</PLAQUETTETERM>
```

- The term is any operator expression over the site basis operators and the
  bond/site operators of the library, like a `BONDTERM`. The number of roles
  must match the plaquette size. Plaquettes may have 2, 3, 4 or more sites.
- Without a `type` attribute the term applies to every plaquette type, and `#`
  is replaced by the type, as for site and bond terms. A plaquette type with no
  matching term contributes nothing.
- Fermion signs follow the operator order in the term string and the ALPS
  site order (Jordan-Wigner over sites in index order, up before down on a
  site), the same convention as bond terms.
- Terms need not conserve particle number. If a term maps a basis state out of
  the sector set by `CONSERVED_QUANTUMNUMBERS`, the run stops with an error
  instead of dropping the matrix element.

A worked example (one 2x2 plaquette followed through the local table and the
matrix fill), the benchmarks against QDK CASCI and exact results, and the
scripts that produced them are in [`n_body_terms/`](n_body_terms/README.md).
