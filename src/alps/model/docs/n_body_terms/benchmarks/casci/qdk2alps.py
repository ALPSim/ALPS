"""QDK Hamiltonian -> ALPS input: a custom lattice and a Hamiltonian former.

The reverse of alps2qdk.  Reads the integrals back out of a QDK Hamiltonian and
writes what ALPS needs to diagonalise it:

    Hamiltonian -> <GRAPH> lattice XML
                -> <HAMILTONIAN> model XML   (only when ALPS has no such model)
                -> parameter file
                -> fulldiag / sparsediag

Two structural facts drive the whole design:

* **ALPS parameters are per type, not per site or bond.**  A model says
  ``-t#*fermion_hop(i,j)`` and ALPS resolves ``t#`` to ``t0``/``t1``/... from
  the *bond type*.  So every distinct hopping value needs its own edge type and
  every distinct on-site energy its own vertex type.  That is what
  ``assign_types`` does, and it is why an arbitrary dense h_pq degenerates into
  one type per matrix element.
* **ALPS graphs are 1-indexed.**  QDK sites are 0-indexed.  Everything written
  here shifts by one; ``site_of`` / ``vertex_of`` are the only places that
  conversion happens.

Scope: lattice-shaped Hamiltonians -- Hubbard, Huckel, PPP -- go through
SITETERM/BONDTERM as before.  Anything else (a molecular Hamiltonian's dense
(pq|rs), pair hopping, correlated hopping, 3- and 4-orbital elements) goes
through ALPS's n-site terms: every set of >= 2 distinct orbitals touched by a
nonzero (pq|rs) becomes an <NSITE> of the graph with its own type, and the
model gets one <NSITETERM> per type holding those elements with numeric
coefficients.  Nothing is truncated.  Needs an ALPS build with NSITETERM.
"""

from __future__ import annotations

import numpy as np

__all__ = ["assign_types", "lattice_xml", "model_xml", "general_model_xml",
           "nsite_groups", "parm_text", "to_alps"]

TOL = 1e-10

# former_for returns this name when the two-body part needs n-site terms.
GENERAL_MODEL = "qdk general fermion"


# ---------------------------------------------------------------- type folding

def assign_types(values, tol: float = TOL):
    """Fold near-equal values onto ALPS type indices.

    Values may be scalars or tuples.  Tuples matter: a bond type carries *all*
    of its couplings at once (t#, V#), so two bonds share a type only when
    every coupling agrees.  Clustering t alone would silently merge bonds that
    differ in V.

    Returns (type_of_each_value, distinct_values), ordered by first appearance,
    so a uniform lattice gets the single type 0.
    """
    def same(a, b):
        if isinstance(a, tuple):
            return all(abs(x - y) <= tol for x, y in zip(a, b))
        return abs(a - b) <= tol

    distinct: list = []
    types: list[int] = []
    for v in values:
        for k, d in enumerate(distinct):
            if same(v, d):
                types.append(k)
                break
        else:
            types.append(len(distinct))
            distinct.append(v)
    return types, distinct


def edge_couplings(ham, h1: np.ndarray, tol: float = TOL):
    """Per-edge (p, q, t, V) read off the integrals.

    t = -h_pq        (ALPS writes hopping as -t)
    V = (pp|qq)      density-density

    Exchange (pq|qp) is not read here; former_for sends it to n-site groups.
    """
    n = h1.shape[0]
    g = None
    if ham.has_two_body_integrals():
        g = np.asarray(ham.get_two_body_integrals()[0]).reshape(n, n, n, n)

    out = []
    for p in range(n):
        for q in range(p + 1, n):
            t = -float(h1[p, q])
            V = float(g[p, p, q, q]) if g is not None else 0.0
            if abs(t) > tol or abs(V) > tol:
                out.append((p, q, t, V))
    return out


def site_couplings(ham, h1: np.ndarray):
    """Per-site (mu, U).

    mu = -h_pp     (ALPS writes the site energy as -mu)
    U  = (pp|pp)   on-site Coulomb
    """
    n = h1.shape[0]
    u = _onsite_u(ham, n)
    if u is None:
        u = [0.0] * n
    return [(-float(h1[p, p]), float(u[p])) for p in range(n)]


# ---------------------------------------------------------------- lattice XML

def lattice_xml(ham, h1: np.ndarray, name: str = "qdk-graph", tol: float = TOL,
                groups: dict | None = None) -> tuple[str, dict]:
    """Integrals -> an ALPS <GRAPH>, plus the parameters its types imply.

    A bond carries (t, V) and is typed on both together; a site carries
    (mu, U) and is likewise typed on both.  Typing a site on mu alone is a bug:
    a lattice with uniform on-site energy but per-site U then emits no <VERTEX>
    elements, every site defaults to type 0, ALPS resolves U# to U0 everywhere
    and silently solves a uniform-U model.
    """
    h1 = np.asarray(h1, dtype=float)
    n = h1.shape[0]

    edges = edge_couplings(ham, h1, tol)
    if groups is not None:
        # General mode: every two-body element over >= 2 sites lives in a
        # n-site group, so the edges carry the hopping only.
        edges = [(p, q, t, 0.0) for p, q, t, _ in edges if abs(t) > tol]
    etypes, evals = assign_types([(t, V) for _, _, t, V in edges], tol)
    vtypes, vvals = assign_types(site_couplings(ham, h1), tol)

    lines = [f'<GRAPH name="{name}" vertices="{n}">']
    # Emit VERTEX elements only when the sites are not all alike; ALPS defaults
    # every vertex to type 0, so a uniform lattice needs none.
    if len(vvals) > 1:
        for p in range(n):
            lines.append(f'  <VERTEX id="{p + 1}" type="{vtypes[p]}"/>')
    for (p, q, _, _), et in zip(edges, etypes):
        lines.append(f'  <EDGE type="{et}" source="{p + 1}" target="{q + 1}"/>')
    # One n-site group type per orbital set: its coefficients are its own.
    for k, sites in enumerate(groups or {}):
        lines.append(f'  <NSITE type="{k}">'
                     + "".join(f'<SITE id="{s + 1}"/>' for s in sites) + "</NSITE>")
    lines.append("</GRAPH>")

    params: dict = {}
    for k, (t, V) in enumerate(evals):
        params[f"t{k}"] = t
        if abs(V) > tol:
            params[f"V{k}"] = V
    if len(vvals) > 1:
        for k, (mu, U) in enumerate(vvals):
            if abs(mu) > tol:
                params[f"mu{k}"] = mu
            if abs(U) > tol:
                params[f"U{k}"] = U
    elif vvals:
        mu, U = vvals[0]
        if abs(mu) > tol:
            params["mu"] = mu
        if abs(U) > tol:
            params["U"] = U
    return "\n".join(lines), params


# ------------------------------------------------------------------ model XML

def _onsite_u(ham, n: int, tol: float = TOL):
    """(pp|pp) per site, or None when there is no two-body part."""
    if not ham.has_two_body_integrals():
        return None
    g = np.asarray(ham.get_two_body_integrals()[0]).reshape(n, n, n, n)
    return [float(g[p, p, p, p]) for p in range(n)]


def former_for(ham, h1: np.ndarray):
    """Decide which ALPS model expresses this Hamiltonian.

    Returns (model_name, needs_custom_xml).  Raises when the integrals have no
    lattice-model form at all.
    """
    n = h1.shape[0]
    u = _onsite_u(ham, n)

    # A one-body-only Hamiltonian is still SPINFUL: the restricted container
    # shares h_pq between alpha and beta. ALPS's `spinless fermions` has no
    # Nup/Ndown quantum numbers, so it would ignore the (n_alpha, n_beta)
    # sector and return one spin species' energy -- exactly half the right
    # answer. Use the spinful model with U left at its default of 0.
    if not ham.has_two_body_integrals() or u is None:
        return "fermion Hubbard", False

    g = np.asarray(ham.get_two_body_integrals()[0]).reshape(n, n, n, n)

    # Only skip the scan when the two-body block is genuinely empty. Testing the
    # on-site U alone is not enough: a Hamiltonian with zero (pp|pp) can still
    # carry 4-index terms, and returning here would accept them unexamined and
    # then drop them.
    if not np.any(np.abs(g) > TOL):
        return "fermion Hubbard", False

    # What ALPS can name is fixed by term arity: SITETERM is 1 site, BONDTERM
    # is 2.  So the only question about a (pq|rs) element is how many *distinct*
    # sites it touches -- not which pattern it matches.  Density-density
    # (p,p,q,q) touches 2 and is a BONDTERM; anything else goes into
    # n-site groups.
    idx = np.array(np.nonzero(np.abs(g) > TOL)).T
    needs_nsite = False
    for row in idx:
        p, q, r, s = (int(x) for x in row)
        sites = {p, q, r, s}

        # Arity: 3 or 4 distinct sites need ALPS's NSITETERM.
        if len(sites) > 2:
            needs_nsite = True
            continue
        if len(sites) == 1:
            continue                                   # on-site U, read by _onsite_u

        # Two sites is necessary but not sufficient: the bond model has only
        # density-density (pp|qq). Exchange (pq|qp) is not ALPS's
        # exchange(i,j), which is S_i.S_j and lacks the density part, so it
        # goes into 2-site groups with pair hopping (pq|pq) and correlated
        # hopping (pp|pq).
        if not (p == q and r == s):
            needs_nsite = True

    if needs_nsite:
        return GENERAL_MODEL, True

    _, uvals = assign_types(u)
    # ALPS's stock `fermion Hubbard` covers a single uniform U and nothing else.
    if len(uvals) > 1:
        return "qdk fermion Hubbard", True
    return "fermion Hubbard", False


def model_xml(name: str = "qdk fermion Hubbard") -> str:
    """A Hamiltonian former, for when ALPS's own library has no match.

    Mirrors ALPS's `fermion Hubbard` but leaves every coefficient as a '#'
    parameter so it resolves per site/bond type.
    """
    return f"""<HAMILTONIAN name="{name}">
  <PARAMETER name="mu" default="0"/>
  <PARAMETER name="t" default="1"/>
  <PARAMETER name="U" default="0"/>
  <PARAMETER name="V" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">
    <PARAMETER name="mu#" default="mu"/>
    <PARAMETER name="U#" default="U"/>
    -mu#*n(i)+U#*n_up(i)*n_down(i)
  </SITETERM>
  <BONDTERM source="i" target="j">
    <PARAMETER name="t#" default="0"/>
    <PARAMETER name="V#" default="0"/>
    -t#*fermion_hop(i,j) + V#*n(i)*n(j)
  </BONDTERM>
</HAMILTONIAN>"""


def nsite_groups(ham, n: int, tol: float = TOL) -> dict:
    """Two-body elements grouped by the orbitals they touch.

    Returns {sorted tuple of >= 2 distinct sites: [(p, q, r, s, (pq|rs)), ...]}
    over every nonzero element; on-site (pp|pp) is left to the site term.
    """
    groups: dict = {}
    if not ham.has_two_body_integrals():
        return groups
    g = np.asarray(ham.get_two_body_integrals()[0]).reshape(n, n, n, n)
    for p, q, r, s in np.array(np.nonzero(np.abs(g) > tol)).T:
        p, q, r, s = int(p), int(q), int(r), int(s)
        key = tuple(sorted({p, q, r, s}))
        if len(key) > 1:
            groups.setdefault(key, []).append((p, q, r, s, float(g[p, q, r, s])))
    return groups


def _nsite_term(sites: tuple, elements: list) -> str:
    """The elements of one orbital set as an ALPS term string.

    Chemist's notation, as everywhere in this module:
        H2 = 1/2 sum (pq|rs) c+_{p,sigma} c+_{r,tau} c_{s,tau} c_{q,sigma}.
    Roles r0, r1, ... name the n-site group's sites in order.
    """
    role = {s: f"r{k}" for k, s in enumerate(sites)}
    out = []
    for p, q, r, s, v in elements:
        c = 0.5 * v
        for a in ("up", "down"):
            for b in ("up", "down"):
                if p == r and a == b:
                    continue                    # cdag cdag on one mode: zero
                if q == s and a == b:
                    continue
                out.append(f"({c!r})*cdag_{a}({role[p]})*cdag_{b}({role[r]})"
                           f"*c_{b}({role[s]})*c_{a}({role[q]})")
    return "\n    + ".join(out) if out else "0"


def general_model_xml(groups: dict, name: str = GENERAL_MODEL) -> str:
    """Hamiltonian former with SITETERM, BONDTERM and one NSITETERM per
    n-site group type (= per orbital set in ``groups``, in the same order as
    lattice_xml emits the n-site groups)."""
    terms = []
    for k, (sites, elements) in enumerate(groups.items()):
        roles = " ".join(f"r{j}" for j in range(len(sites)))
        terms.append(f'  <NSITETERM type="{k}" sites="{roles}">\n    '
                     + _nsite_term(sites, elements) + "\n  </NSITETERM>")
    return f"""<HAMILTONIAN name="{name}">
  <PARAMETER name="mu" default="0"/>
  <PARAMETER name="t" default="1"/>
  <PARAMETER name="U" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">
    <PARAMETER name="mu#" default="mu"/>
    <PARAMETER name="U#" default="U"/>
    -mu#*n(i)+U#*n_up(i)*n_down(i)
  </SITETERM>
  <BONDTERM source="i" target="j">
    <PARAMETER name="t#" default="0"/>
    -t#*fermion_hop(i,j)
  </BONDTERM>
""" + "\n".join(terms) + "\n</HAMILTONIAN>"


# ------------------------------------------------------------- parameter file

def parm_text(params: dict) -> str:
    def fmt(v):
        return f'"{v}"' if isinstance(v, str) else repr(v)
    return "\n".join(f"{k}={fmt(v)}" for k, v in params.items()) + "\n"


def to_alps(ham, outdir: str, name: str = "qdk-graph", extra: dict | None = None) -> dict:
    """Write lattice XML, model XML (if needed) and a parameter file.

    Returns the paths written plus the parameter dict, so a caller can hand
    them straight to fulldiag.
    """
    import os

    # Only the alpha block is read below, which is correct exactly when the
    # Hamiltonian is restricted (alpha and beta share h_pq). An unrestricted
    # one would silently lose its beta block, so refuse instead.
    if ham.has_orbitals() and not ham.get_orbitals().is_restricted():
        raise NotImplementedError(
            "unrestricted Hamiltonian: qdk2alps reads only the alpha one-body "
            "block, and ALPS's spinful fermion models share one set of "
            "couplings between the spin species. Emitting this would silently "
            "drop the beta integrals."
        )

    h1_a, h1_b = ham.get_one_body_integrals()
    h1 = np.asarray(h1_a, dtype=float)
    if not np.allclose(h1, np.asarray(h1_b, dtype=float), atol=1e-10):
        raise NotImplementedError(
            "alpha and beta one-body integrals differ; ALPS's fermion Hubbard "
            "applies one t/mu to both spin species."
        )
    n = h1.shape[0]

    model, custom = former_for(ham, h1)
    groups = nsite_groups(ham, n) if model == GENERAL_MODEL else None
    graph, params = lattice_xml(ham, h1, name, groups=groups)

    # mu and U are emitted by lattice_xml, typed jointly per site.
    params["LATTICE"] = name
    params["MODEL"] = model
    if extra:
        params.update(extra)

    os.makedirs(outdir, exist_ok=True)
    out = {}

    lat_path = os.path.join(outdir, "lattices.xml")
    with open(lat_path, "w") as fh:
        fh.write("<LATTICES>\n" + graph + "\n</LATTICES>\n")
    out["lattices"] = lat_path
    params["LATTICE_LIBRARY"] = lat_path

    if custom:
        # A MODEL_LIBRARY replaces models.xml wholesale, so the emitted model's
        # <BASIS ref="fermion"/> would dangle. Carry the stock library's body
        # across and append the new Hamiltonian to it.
        # the repository's stock library; this file is
        # src/alps/model/docs/n_body_terms/benchmarks/casci/qdk2alps.py
        base = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            *[".."] * 7, "src", "alps", "resources", "models.xml")
        with open(base) as fh:
            body = fh.read()
        body = body[body.index("<MODELS>") + len("<MODELS>"):body.rindex("</MODELS>")]

        mod_path = os.path.join(outdir, "models.xml")
        with open(mod_path, "w") as fh:
            former = (general_model_xml(groups, model) if model == GENERAL_MODEL
                      else model_xml(model))
            fh.write("<MODELS>\n" + body + "\n" + former + "\n</MODELS>\n")
        out["models"] = mod_path
        params["MODEL_LIBRARY"] = mod_path

    parm_path = os.path.join(outdir, "parm")
    with open(parm_path, "w") as fh:
        fh.write(parm_text(params))
    out["parm"] = parm_path
    out["params"] = params
    return out
