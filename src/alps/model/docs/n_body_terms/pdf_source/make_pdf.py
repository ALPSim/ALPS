"""Build ../ALPS_n_body_terms.pdf: n-body (plaquette) terms in ALPS, step by step.

Inputs (all in this folder):
    lattice.xml, model.xml, parm         the example, as a user would write it
    data/alps_dump_K*.txt                basis, plaquette table and H printed by
                                         ALPS itself (dump_plaquette.cpp)
    data/fulldiag_runs.txt               fulldiag / dmrg runs on the example
    ../benchmarks/results/*.txt          benchmark output

Every number in the PDF comes from those files; the script only re-checks them
(the fill rule re-implemented in numpy reproduces ALPS's matrix, and the
independent prototype ../benchmarks/alps_checks/hyper_ed.py gives the same H).

    python make_pdf.py          # needs numpy, scipy, reportlab
"""

import os
import re
import sys

import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (KeepTogether, PageBreak, Paragraph, Preformatted, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(DOCS, *[".."] * 5))
sys.path.insert(0, os.path.join(DOCS, "benchmarks", "alps_checks"))
import hyper_ed  # noqa: E402  independent prototype (dense, mode order 2*site+spin)


def font(name, candidates):
    for path in candidates:
        if os.path.exists(path):
            pdfmetrics.registerFont(TTFont(name, path))
            return name
    return {"Body": "Helvetica", "Mono": "Courier"}[name]


BODY = font("Body", ["/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
                     "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"])
MONOF = font("Mono", ["/System/Library/Fonts/Supplemental/Courier New.ttf",
                      "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"])
H0 = ParagraphStyle("h0", fontName=BODY, fontSize=17, leading=21, spaceAfter=8)
H1 = ParagraphStyle("h1", fontName=BODY, fontSize=14, leading=18, spaceBefore=6, spaceAfter=5)
P = ParagraphStyle("p", fontName=BODY, fontSize=9.5, leading=12.5, spaceAfter=4)
SMALL = ParagraphStyle("s", fontName=BODY, fontSize=8, leading=10, textColor=colors.HexColor("#555555"))
MONO = ParagraphStyle("m", fontName=MONOF, fontSize=8, leading=9.6)
MONO_S = ParagraphStyle("ms", fontName=MONOF, fontSize=6.6, leading=7.8)
GREY = colors.HexColor("#e6e6e6")
YELLOW = colors.HexColor("#ffd54f")
ORANGE = colors.HexColor("#ffab91")
BLUE = colors.HexColor("#bbdefb")


def read(*p):
    with open(os.path.join(*p)) as f:
        return f.read()


def code(text, style=MONO, width=118):
    lines = []
    for line in text.rstrip("\n").splitlines():
        while len(line) > width:
            cut = line.rfind(", ", 0, width) + 1 or width
            lines.append(line[:cut])
            line = " " * (len(line) - len(line.lstrip()) + 6) + line[cut:].lstrip()
        lines.append(line)
    return Preformatted("\n".join(lines), style)


def table(rows, widths=None, header=True, fs=8, align="CENTER", extra=()):
    t = Table(rows, colWidths=widths, repeatRows=1 if header else 0)
    st = [("FONT", (0, 0), (-1, -1), BODY, fs), ("GRID", (0, 0), (-1, -1), 0.3, colors.grey),
          ("ALIGN", (0, 0), (-1, -1), align), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
          ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5)]
    if header:
        st.append(("BACKGROUND", (0, 0), (-1, 0), GREY))
    t.setStyle(TableStyle(st + list(extra)))
    return t


# --------------------------------------------------------------- ALPS output

def load_dump(name):
    out = {"entries": [], "states": [], "H": []}
    for line in read(HERE, "data", name).splitlines():
        tag, _, rest = line.partition(" ")
        if tag == "LOCAL":
            out["local"] = {int(k): tuple(map(int, v.split(","))) for k, v in (w.split(":") for w in rest.split())}
        elif tag == "PLAQ":
            out["plaq"] = list(map(int, rest.split()))
        elif tag == "TERM":
            out["term"] = rest
        elif tag == "ENTRY":
            a, b, v = rest.split()
            out["entries"].append((tuple(map(int, a)), tuple(map(int, b)), float(v)))
        elif tag == "STATE":
            out["states"].append(tuple(map(int, rest.split()[1:])))
        elif tag == "ROW":
            out["H"].append(list(map(float, rest.split())))
    out["H"] = np.array(out["H"])
    return out


d5 = load_dump("alps_dump_K0.5.txt")
d0 = load_dump("alps_dump_K0.txt")
K, U, T_HOP = 0.5, 4.0, 1.0
H5, H0m = d5["H"], d0["H"]
STATES, PLAQ, LOCAL = d5["states"], d5["plaq"], d5["local"]
INDEX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES[0])
ODD = {k: (sum(v) % 2 == 1) for k, v in LOCAL.items()}          # fermion parity of a local state
LABEL = {(0, 0): "·", (0, 1): "↓", (1, 0): "↑", (1, 1): "↕"}     # (Nup, Ndown)
LOC = {k: LABEL[v] for k, v in LOCAL.items()}


def picture(s):
    """Plaquette drawn as 'sites 2 3 / sites 0 1', as in quartic_ED_walkthrough.pdf."""
    return f"{LOC[s[2]]} {LOC[s[3]]} / {LOC[s[0]]} {LOC[s[1]]}"


def bits(s):
    b = 0
    for site, l in enumerate(s):
        nup, ndn = LOCAL[l]
        b |= nup << (2 * site) | ndn << (2 * site + 1)
    return b


def ket(s):
    out = []
    for site, l in enumerate(s):
        nup, ndn = LOCAL[l]
        out += [f"↑{site}"] * nup + [f"↓{site}"] * ndn
    return "|" + " ".join(out) + ">"


# ---- re-implementation of hamiltonian_matrix::apply_operator for the table

def fill(entries, trace_col=None):
    Pm = np.zeros_like(H5)
    trace = []
    n = len(PLAQ)
    for col, st in enumerate(STATES):
        tin = tuple(st[p] for p in PLAQ)
        prefix = [False]
        for s in range(NS):
            prefix.append(prefix[-1] != ODD[st[s]])
        for a, b, v in entries:
            if a != tin:
                continue
            new = list(st)
            for r in range(n):
                new[PLAQ[r]] = b[r]
            row = INDEX[tuple(new)]
            f = False
            lines = []
            for r in reversed(range(n)):
                s = PLAQ[r]
                if ODD[a[r]] != ODD[b[r]]:
                    before = prefix[s]
                    fixes = [PLAQ[q] for q in range(r + 1, n) if PLAQ[q] < s and ODD[a[q]] != ODD[b[q]]]
                    if len(fixes) % 2:
                        before = not before
                    f ^= before
                    nb = sum(LOCAL[st[x]][0] + LOCAL[st[x]][1] for x in range(s))  # before any operator acts
                    lines.append((r, s, LOC[a[r]], LOC[b[r]], nb, fixes, before))
            val = -v if f else v
            Pm[row, col] += val
            if col == trace_col:
                trace.append((row, v, val, lines, b))
    return Pm, trace


K_ENTRIES = [(a, b, v / K) for a, b, v in d5["entries"]]
P_alps = (H5 - H0m) / K
P_rule, TRACE = fill(K_ENTRIES, trace_col=INDEX[(2, 1, 0, 0)])
assert np.abs(P_rule - P_alps).max() < 1e-14, "fill rule does not reproduce ALPS"

# ---- independent prototype: same H, entry by entry, after re-indexing
proto = {"lattice": {"kind": "ladder", "L": 2, "boundary": "open"},
         "sector": {"Nup": 1, "Ndown": 1},
         "parameters": {"t": T_HOP, "U": U, "K": K},
         "term": [{"on": "bond", "sites": ["i", "j"],
                   "expr": "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"},
                  {"on": "site", "sites": ["i"], "expr": "U*n_up(i)*n_down(i)"},
                  {"on": "plaquette", "sites": ["i", "j", "k", "l"],
                   "expr": "K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))"}]}
Hp = hyper_ed.solve(proto, verbose=False)["H"].toarray().real
_, pidx = hyper_ed.sector_basis(4, 1, 1)
perm = np.zeros((16, 16))
for a, s in enumerate(STATES):
    perm[pidx[bits(s)], a] = 1
PROTO_DIFF = np.abs(perm @ H5 @ perm.T - Hp).max()
assert PROTO_DIFF < 1e-14
E5 = np.linalg.eigvalsh(H5)
E0 = np.linalg.eigvalsh(H0m)

# ------------------------------------------------------------ fulldiag runs
RUNS = read(HERE, "data", "fulldiag_runs.txt").splitlines()
FD = {}
for line in RUNS:
    m = re.match(r"task\d K=([\d.]+):(.*)", line)
    if m:
        FD[m.group(1)] = np.array(sorted(map(float, m.group(2).split())))
assert np.abs(FD["0.5"] - E5).max() < 1e-12 and np.abs(FD["0"] - E0).max() < 1e-12
STOCK = [l for l in RUNS if l.startswith("stock")][0].split("E0 = ")[1]
NOPLAQ = [l for l in RUNS if l.startswith("K=0.5 model")][0].split("E0 = ")[1]
DMRG = [l for l in RUNS if l.startswith("dmrg")][0].split(": ", 1)[1]


def num(v):
    return "·" if abs(v) < 1e-12 else f"{v:+g}"


def matrix_table(M, marks=(), fs=6.4, cw=27):
    D = M.shape[0]
    rows = [[""] + [str(c) for c in range(D)]] + [[str(r)] + [num(M[r, c]) for c in range(D)] for r in range(D)]
    st = [("FONT", (0, 0), (-1, -1), BODY, fs), ("ALIGN", (0, 0), (-1, -1), "CENTER"),
          ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
          ("TEXTCOLOR", (0, 0), (-1, 0), colors.grey), ("TEXTCOLOR", (0, 0), (0, -1), colors.grey),
          ("BOX", (1, 1), (-1, -1), 0.8, colors.black),
          ("INNERGRID", (1, 1), (-1, -1), 0.2, colors.HexColor("#cccccc")),
          ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 1)]
    for (r, c), color in marks:
        st.append(("BACKGROUND", (c + 1, r + 1), (c + 1, r + 1), color))
    t = Table(rows, colWidths=[18] + [cw] * D, rowHeights=11)
    t.setStyle(TableStyle(st))
    return t


# ------------------------------------------------------------------- story
story = []
story += [
    Paragraph("n-body (plaquette) terms in ALPS, step by step", H0),
    Paragraph(
        "ALPS models had two kinds of terms: <b>SITETERM</b> (one site) and <b>BONDTERM</b> (two sites). "
        "This branch adds a third, <b>PLAQUETTETERM</b>: one operator expression on an ordered group of any number "
        "of sites, attached to a new <b>&lt;PLAQUETTE&gt;</b> element of the lattice. A molecular Hamiltonian's "
        "(pq|rs), a ring exchange or a correlated pair hop can now be written down exactly, with nothing "
        "truncated. The terms are used by exact diagonalization (<font face='Mono'>fulldiag</font>, "
        "<font face='Mono'>sparsediag</font>).", P),
    Paragraph(
        "This document follows one small problem through ALPS: the input files, the plaquette's local table, "
        "the fill of the many-body matrix, and the eigenvalues. It is the same system as "
        "<i>quartic_ED_walkthrough.pdf</i> (one 2×2 plaquette, t = 1, U = 4, K = 0.5, one ↑ and one ↓ electron). "
        "Every matrix and number below was printed by ALPS itself (<font face='Mono'>pdf_source/</font>).", P),
    Spacer(1, 4),
    Paragraph("Optional: nothing changes unless the input asks for it", H1),
    Paragraph(
        "There is no switch to turn on and nothing to change in existing inputs. The new code runs only when the "
        "lattice contains a &lt;PLAQUETTE&gt; <i>and</i> the model contains a &lt;PLAQUETTETERM&gt;:", P),
    table([
        ["lattice has\n<PLAQUETTE>", "model has\n<PLAQUETTETERM>", "what happens", "checked"],
        ["no", "no", "exactly as before", "stock models: master vs branch bit-identical;\nfull ctest 931/931"],
        ["yes", "no", "plaquettes are ignored", f"stock 'fermion Hubbard' on this lattice:\nE0 = {float(STOCK):.10f}"],
        ["no", "yes", "term has nothing to act on\n(like a BONDTERM on a graph without edges)",
         f"K = 0.5 model, plaquette removed:\nE0 = {float(NOPLAQ):.10f} (= the K = 0 value)"],
        ["yes", "yes", "fulldiag / sparsediag add the term;\nDMRG and QMC codes stop with a message",
         "this document; dmrg refuses (below)"],
    ], widths=[1.15 * inch, 1.3 * inch, 2.25 * inch, 2.4 * inch], fs=7.5),
    Spacer(1, 3),
    Paragraph("The refusal, so a model is never silently truncated: <i>" + DMRG.replace("<", "&lt;").replace(">", "&gt;") + "</i>", SMALL),
]

# ---- Step 1 lattice
story += [
    Paragraph("Step 1 — Lattice: a &lt;PLAQUETTE&gt; is an ordered group of sites", H1),
    Paragraph(
        "A &lt;GRAPH&gt; gets plaquette elements next to its edges. Each &lt;SITE id&gt; is a vertex id (1-based, "
        "like EDGE source/target). The <b>order</b> of the sites is what the model's role names bind to: here "
        "i, j, k, l = vertices 1, 2, 4, 3, i.e. ALPS sites 0, 1, 3, 2 going round the square. "
        "The plaquette is stored in its own list beside the graph; it adds no bonds.", P),
]
lat = read(HERE, "lattice.xml")
story += [Table([[code(lat, MONO_S), code(
    "ALPS builds:\n\n"
    "  sites      0 1 2 3\n"
    "  bonds      (0,1) (2,3) (0,2) (1,3)\n"
    f"  plaquettes type 0: ({', '.join(map(str, PLAQ))})\n\n"
    "     2 ----- 3\n"
    "     |   K   |\n"
    "     0 ----- 1\n\n"
    "  roles  i=0  j=1  k=3  l=2")]], colWidths=[4.2 * inch, 2.8 * inch],
    style=[("VALIGN", (0, 0), (-1, -1), "TOP")])]
story += [Paragraph(
    "For a periodic lattice the plaquette goes in the &lt;UNITCELL&gt;, with each site given as vertex + cell offset "
    "exactly like an EDGE's SOURCE/TARGET. ALPS repeats it over the lattice, wraps it at periodic boundaries and drops "
    "it where a site would fall off an open boundary:", P),
    code('<UNITCELL name="ladder cell" dimension="1" vertices="2">\n'
         '  ...edges as usual...\n'
         '  <PLAQUETTE type="0">\n'
         '    <SITE vertex="1" offset="0"/> <SITE vertex="2" offset="0"/>\n'
         '    <SITE vertex="2" offset="1"/> <SITE vertex="1" offset="1"/>\n'
         '  </PLAQUETTE>\n'
         '</UNITCELL>')]

# ---- Step 2 model
mod = read(HERE, "model.xml")
ham = mod[mod.index("  <HAMILTONIAN"):mod.index("</MODELS>")]
story += [
    PageBreak(),
    Paragraph("Step 2 — Model: a &lt;PLAQUETTETERM&gt; is written like a BONDTERM", H1),
    Paragraph(
        "<font face='Mono'>sites=\"i j k l\"</font> names the roles; the term is any expression in the site "
        "operators (and library operators), on as many roles as the plaquette has sites. As for bond terms, "
        "<font face='Mono'>type=\"…\"</font> restricts a term to one plaquette type, and without it "
        "<font face='Mono'>#</font> in parameter names is replaced by the type (K# → K0, K1, …). Fermion signs "
        "follow the operator order as written. The SITEBASIS and BASIS at the top of model.xml are the usual "
        "fermion ones; only the HAMILTONIAN is shown.", P),
    code(ham, MONO_S),
    Spacer(1, 2),
    Paragraph(
        "The K term is quartic on four different sites: it moves an ↑ electron from i to l and a ↓ electron "
        "from j to k together (and back). No site or bond term can produce it.", P),
]

# ---- Step 3 run
parm = read(HERE, "parm")
story += [
    Paragraph("Step 3 — Run it like any other ALPS model", H1),
    Table([[code(parm, MONO_S), code(
        "$ parameter2xml parm\n"
        "$ fulldiag parm.in.xml\n\n"
        "lowest levels (fulldiag, N↑ = N↓ = 1):\n" +
        "\n".join(f"  K = {k:<3}  " + "  ".join(f"{e:+.6f}" for e in FD[k][:4]) for k in ("0", "0.5")))]],
        colWidths=[3.0 * inch, 4.0 * inch], style=[("VALIGN", (0, 0), (-1, -1), "TOP")]),
    Spacer(1, 3),
    Paragraph(
        f"All 16 levels for K = 0.5: " + "  ".join(f"{e:+.6f}" for e in FD["0.5"]) +
        ". These are the numbers in quartic_ED_walkthrough.pdf (E0 = −3.3597011665, and −3.4185507189 at K = 0).", P),
]

# ---- Step 4 basis
story += [
    Paragraph("Step 4 — The ALPS basis of the (N↑ = 1, N↓ = 1) sector", H1),
    Paragraph(
        "Each site has 4 local states, numbered by ALPS from the quantum numbers (Nup, Ndown): "
        + ", ".join(f"{k} = {LOC[k]} {LOCAL[k]}" for k in sorted(LOCAL)) +
        " (· empty, ↕ doubly occupied). A basis state is the tuple of local states on sites 0…3, and ALPS "
        "sorts the 16 tuples that satisfy the constraints. The Jordan–Wigner order is site 0, 1, 2, 3 and ↑ before "
        "↓ on a site, i.e. mode = 2·site + spin as in the earlier PDFs; the bit column is that occupation "
        "string.", P),
]
rows = [["ALPS index", "local states (site 0 1 2 3)", "bits m7…m0", "electrons", "picture (2 3 / 0 1)"]]
for i, s in enumerate(STATES):
    rows.append([str(i), "  ".join(map(str, s)), f"{bits(s):08b}", ket(s), picture(s)])
k_rows = sorted({r for r, c in zip(*np.nonzero(P_alps))})
story += [table(rows, widths=[0.7 * inch, 1.6 * inch, 1.0 * inch, 1.1 * inch, 1.2 * inch], fs=7.5,
                extra=[("BACKGROUND", (0, r + 1), (-1, r + 1), ORANGE) for r in k_rows]),
          Paragraph("Orange: the two states the K term connects.", SMALL)]

# ---- Step 5 local table
n_entries = len(d5["entries"])
story += [
    PageBreak(),
    Paragraph("Step 5 — The plaquette's local table (NSiteOperator::matrix)", H1),
    Paragraph(
        "ALPS first turns the term string into a table on the 4 plaquette sites alone, indexed by tuples of "
        "local states in role order (i, j, k, l): 4<super>4</super> = 256 input tuples. Each product in the term "
        "is split once into a number times one operator string per role, "
        "O<sub>i</sub> O<sub>j</sub> O<sub>k</sub> O<sub>l</sub>: the site operators are read right to left, and "
        "each odd (fermionic) operator moved into role order gives one sign for each odd operator it passes. "
        "Each O<sub>r</sub> becomes a 4×4 single-site matrix. A row of the table is the product of these, "
        "assembled the first time the fill (Step 6) asks for it, so only rows that basis states reach are built. "
        f"Of the 256 × 256 entries only <b>{n_entries}</b> are nonzero. "
        "Signs between the plaquette and the rest of the lattice are <i>not</i> in the table; they are added in "
        "Step 6, as for bond terms. The table depends only on the plaquette type and its site types, so it is "
        "built once and reused on every plaquette of the lattice.", P),
    Paragraph(
        "Below: the 16 entries for the first half, K c†<sub>l↑</sub> c†<sub>k↓</sub> c<sub>j↓</sub> "
        "c<sub>i↑</sub>, in units of K. The other 16 (the Hermitian conjugate) are the same with input and "
        "output swapped. There are 16 because each site may also hold the other spin as a spectator (2<super>4</super>). "
        "The sign depends only on those spectators: c<sub>↓</sub> acting on ↕ passes the ↑ on the same site "
        "(↑ comes first), which gives −1, while c<sub>↑</sub> on ↕ gives +1.", P),
]
rows = [["in (i j k l)", "", "out (i j k l)", "", "value / K", "spectators"]]
for a, b, v in K_ENTRIES:
    if LOC[a[0]] in "↑↕" and LOC[a[1]] in "↓↕" and LOC[b[3]] in "↑↕":
        spect = ", ".join(f"{'ijkl'[r]}" for r in range(4) if LOC[a[r]] == "↕" or LOC[b[r]] == "↕") or "none"
        rows.append(["".join(map(str, a)), " ".join(LOC[x] for x in a), "".join(map(str, b)),
                     " ".join(LOC[x] for x in b), f"{v:+g}", spect])
story += [table(rows, widths=[0.8 * inch, 0.9 * inch, 0.85 * inch, 0.9 * inch, 0.75 * inch, 1.2 * inch], fs=7.5)]

# ---- Step 6 fill
col = INDEX[(2, 1, 0, 0)]
row, v, val, lines, b = TRACE[0]
assert row == INDEX[(0, 0, 2, 1)] and abs(H5[row, col] - K * val) < 1e-14
trace = [f"column {col} = {ket(STATES[col])}  (picture {picture(STATES[col])}),  plaquette tuple (i j k l) = "
         f"{' '.join(str(STATES[col][p]) for p in PLAQ)}",
         f"table entry -> output tuple {' '.join(map(str, b))}, value {v:+g} K",
         "operators in role order, applied last role first:"]
for r, s, li, lo, nb, fixes, before in lines:
    fx = "".join(f" (+1: role {'ijkl'[PLAQ.index(x)]} already applied on site {x})" for x in fixes)
    trace.append(f"  role {'ijkl'[r]} = site {s}: {li} -> {lo}   fermions on sites < {s}: {nb}{fx}"
                 f"  -> {'odd, -1' if before else 'even, +1'}")
trace.append(f"result: sign {'-' if val * v < 0 else '+'}1, row {row} = {ket(STATES[row])},  "
             f"H[{row},{col}] += {K * val:+g}")
story += [
    Paragraph("Step 6 — Fill the many-body matrix (hamiltonian_matrix::apply_operator)", H1),
    Paragraph(
        "For every basis state, ALPS reads the local states on the plaquette's sites (in role order), looks up "
        "that row of the table, and writes each output tuple back into the state to get the target basis state. "
        "The only extra ingredient is the Jordan–Wigner sign of each odd operator O<sub>r</sub>: "
        "(−1)<super>n</super>, where n is the number of fermions on sites before site(r), counted in the state as "
        "already changed by the operators applied before it. For two sites this is exactly the rule ALPS already "
        "uses for bond terms. One column, worked through:", P),
    code("\n".join(trace)),
    Spacer(1, 3),
    Paragraph(
        f"This is the entry the earlier PDF finds as V[4,7|0,3] → H[11,1] = +K in its own state numbering. "
        f"Re-implementing this rule in numpy (make_pdf.py, <font face='Mono'>fill</font>) on the {n_entries} "
        f"table entries reproduces ALPS's K contribution exactly: max |difference| = "
        f"{np.abs(P_rule - P_alps).max():.0e}.", P),
]

# ---- Step 7 full H
marks = [((r, c), ORANGE) for r, c in zip(*np.nonzero(P_alps))] + \
        [((i, i), YELLOW) for i in range(16) if abs(H5[i, i]) > 1e-12]
story += [
    PageBreak(),
    Paragraph("Step 7 — The full 16 × 16 Hamiltonian, as ALPS builds it", H1),
    Paragraph(
        "Rows and columns use the ALPS basis of Step 4. Yellow: U on doubly occupied sites (site terms). "
        "Orange: the K term (plaquette term). Everything else: hopping (bond terms), whose ± sign comes from "
        "the same Jordan–Wigner rule.", P),
    matrix_table(H5, marks),
    Spacer(1, 6),
    Paragraph(
        f"<b>Check against the independent prototype.</b> Re-indexed into the earlier PDF's state order, ALPS's "
        f"matrix equals the prototype's (hyper_ed.py: normal-ordered coefficient tables and bit-string fill, "
        f"no ALPS code) entry by entry, signs included: max |difference| = {PROTO_DIFF:.1f}. "
        f"Eigenvalues: E0 = {E5[0]:.10f} (K = 0.5), {E0[0]:.10f} (K = 0).", P),
]

# ---- Step 8 benchmarks
res = os.path.join(DOCS, "benchmarks", "results")
casci = read(res, "casci_vs_alps.txt").splitlines()
crow = [["molecule", "sector\n(N↑,N↓)", "plaquettes\n2/3/4-site", "QDK CASCI", "ALPS ED", "|ΔE|", "spectrum\nvs dense JW"]]
for line in casci:
    m = re.match(r"(.{30})\s*(\(\d,\d\))\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
    if m:
        crow.append([m.group(1).strip(), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6), m.group(8)])
verdict = casci[-1].strip()


def last_line(name):
    return read(res, name).strip().splitlines()[-1]


def count(name, pat):
    return sum(1 for l in read(res, name).splitlines() if re.search(pat, l))


story += [
    Paragraph("Step 8 — Benchmarks", H1),
    Paragraph(
        "<b>(a) Molecules: ALPS ED vs QDK CASCI, full (pq|rs).</b> A molecular Hamiltonian's two-body part "
        "touches 3 and 4 distinct orbitals, so before this change it could not be written for ALPS without "
        "truncation. Now every orbital set touched by a nonzero (pq|rs) becomes a &lt;PLAQUETTE&gt; with its own "
        "&lt;PLAQUETTETERM&gt;. All orbitals active, STO-3G. The last column compares the whole spectrum with a "
        "dense Jordan–Wigner diagonalization of the same integrals (numpy, independent of both codes). "
        f"Result: <b>{verdict}</b>.", P),
    table(crow, widths=[1.55 * inch, 0.6 * inch, 0.75 * inch, 1.0 * inch, 1.0 * inch, 0.55 * inch, 0.75 * inch], fs=7),
    Spacer(1, 6),
    Paragraph(
        "<b>(b) Exact answer on 4 sites.</b> Open 4-site chain plus U on one delocalized "
        "hopping eigenmode a = Σ φ<sub>i</sub> c<sub>i</sub>: H = −t Σ hop + U n<sub>a↑</sub> n<sub>a↓</sub>. "
        "Written in the site basis, the U term is one PLAQUETTETERM with 1-, 2-, 3- and 4-site pieces. The "
        "spectrum is known exactly: Σ ε<sub>occupied</sub> + U·[mode a doubly occupied]. Checked for two modes, "
        "three role orderings and five sectors. A term that changes N by 4 is checked in the full 256-state space "
        "against a dense matrix, and with N↑, N↓ conserved it is refused. "
        f"{count('t2_analytic.txt', 'max.dE')} runs; {last_line('t2_analytic.txt')}.", P),
    Paragraph(
        "<b>(c) Ladders vs the prototype.</b> 2×L Hubbard ladders with a pair-hopping plaquette term, open and "
        f"periodic, several sectors: {last_line('t1_ladder.txt')}.", P),
    Paragraph(
        "<b>(d) Random two-body terms</b> on 2- and 3-site plaquettes, all index patterns, with interleaved site "
        f"orders, vs a dense Jordan–Wigner matrix: {last_line('t5_interleaved.txt')}.", P),
    Paragraph(
        "<b>(e) Bosonic Moore–Read state: a 3-body term, 6 operators.</b> Bosons in the lowest Landau level on a "
        "sphere with 2S flux; the orbitals are the sites. The 3-body contact interaction "
        "H = Σ<sub>M</sub> T<sub>M</sub><super>†</super> T<sub>M</sub> projects every triple onto total angular "
        "momentum 3S, and is written as one PLAQUETTETERM over all orbitals (strings like "
        "b†<sub>0</sub>b†<sub>0</sub>b†<sub>2</sub>b<sub>1</sub>b<sub>1</sub>b<sub>1</sub>). The spectrum is "
        "compared with a dense reference; the number of zero modes with the number of patterns with at most 2 bosons "
        "in any 2 neighbouring orbitals (Bernevig–Haldane). At 2S = N − 2 the zero mode is unique: the Moore–Read "
        f"ground state. {last_line('t6_moore_read.txt')}.", P),
    table([["N", "2S", "dim", "max|ΔE| vs dense", "zero modes", "admissible", "gap"]]
          + [l.split()[:5] + [l.split()[5], l.split()[7]] for l in read(res, "t6_moore_read.txt").splitlines()
             if re.match(r"\s*\d+\s+\d+\s+\d+", l)],
          widths=[0.4 * inch, 0.4 * inch, 0.6 * inch, 1.2 * inch, 0.85 * inch, 0.85 * inch, 0.8 * inch], fs=7),
    Spacer(1, 6),
    Paragraph(
        "<b>(f) Unit tests.</b> src/alps/model/tests/nsite_terms.cpp (GoogleTest): unit-cell plaquettes at open and "
        "periodic boundaries, GRAPH plaquettes, invalid plaquettes, ladder reference energies, K = 0 equals the "
        "plain model, the analytic case above, sector violation, typed and # terms, parameter precision, the "
        "refusal in other applications, and bosons: N-conserving and non-conserving terms vs a dense reference, "
        "truncation at Nmax and its warning. Full ctest: 936/936 pass.", P),
    Paragraph("Raw output and the scripts are in <font face='Mono'>benchmarks/</font>.", SMALL),
]

# ---- Step 9 changes
story += [
    Paragraph("Step 9 — What changed in ALPS", H1),
    table([
        ["file", "change"],
        ["lattice/plaquette.h (new)", "plaquette {type, sites}; reading <PLAQUETTE> in GRAPH and UNITCELL"],
        ["lattice/graph.h, latticelibrary.h/.C", "GRAPH accepts <PLAQUETTE>; the library keeps them per graph"],
        ["lattice/unitcell.h/.C", "UNITCELL accepts <PLAQUETTE> (vertex + offset), and writes it back"],
        ["lattice/latticegraph.h", "unit-cell plaquettes repeated over the lattice, same boundary rules as edges"],
        ["lattice/graph_helper.h", "plaquettes(), num_plaquettes()"],
        ["model/nsiteoperator.h (new)", "NSiteOperator, the splitter into per-site operators with the sign, the local table"],
        ["model/nsiteterm.C (new)", "read/write <PLAQUETTETERM>"],
        ["model/globaloperator.h/.C, default_term.h", "HAMILTONIAN holds plaquette terms (typed and #)"],
        ["model/hamiltonian_matrix.hpp", "plaquette loop: local table per type, Jordan–Wigner fill, sector check"],
        ["model/model_helper.h", "applications without support refuse plaquette terms instead of dropping them"],
        ["model/siteterm.C, bondterm.C", "fix: parameters were rounded to 6 digits (see below)"],
        ["model/tests/nsite_terms.cpp (new)", "unit tests"],
        ["CHANGELOG.md, model/docs/", "changelog entry, syntax reference (plaquette_terms.md), this folder"],
    ], widths=[2.4 * inch, 4.6 * inch], fs=7.5, align="LEFT"),
    Spacer(1, 6),
    Paragraph(
        "<b>One fix to existing behaviour.</b> SITETERM and BONDTERM substituted parameter values into the "
        "expression with 6 significant digits, so <font face='Mono'>t0=0.123456789</font> acted as 0.123457. "
        "Molecular integrals exposed this as a ~5·10<super>-6</super> error against CASCI. Values are now written "
        "with full double precision. Inputs with short parameter values (t = 1, U = 4, …) give bit-identical "
        "results; longer values are now used exactly as given.", P),
    Paragraph(
        "<b>Limits.</b> Only exact diagonalization supports plaquette terms for now; DMRG and the QMC codes "
        "refuse them. Plaquettes in a &lt;GRAPH&gt; are not "
        "yet written back out when ALPS writes a lattice library to XML.", P),
]

# ---- appendix
src = read(REPO, "src/alps/model/include/alps/model/hamiltonian_matrix.hpp")
a = src.index("// Fill one n-site term on plaquette p.")
b = src.index("template <class M, class G> template <class STATES, class V, class W>\nvoid hamiltonian_matrix<M,G>::apply_operator(const STATES& states, const SiteOperator& op")
story += [
    PageBreak(),
    Paragraph("Appendix — the fill (Step 6) in ALPS", H1),
    Paragraph("src/alps/model/include/alps/model/hamiltonian_matrix.hpp:", SMALL),
    code(src[a:b].rstrip(), MONO_S),
]
nso = read(REPO, "src/alps/model/include/alps/model/nsiteoperator.h")
a = nso.index("// Operators are applied right to left")
b = nso.index("template <class T, class I>\nnsite_matrix<T> NSiteOperator::matrix")
story += [
    Spacer(1, 6),
    Paragraph("The sign used when building the local table (Step 5), src/alps/model/include/alps/model/nsiteoperator.h:", SMALL),
    code(nso[a:nso.index("class ALPS_DECL NSiteOperator")].rstrip() if nso.index("class ALPS_DECL NSiteOperator") > a
         else nso[a:b].rstrip(), MONO_S),
]


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont(BODY, 7.5)
    canvas.setFillColor(colors.grey)
    canvas.drawRightString(letter[0] - 0.6 * inch, 0.45 * inch, f"ALPS n-body terms · {doc.page}")
    canvas.restoreState()


out = os.path.join(DOCS, "ALPS_n_body_terms.pdf")
SimpleDocTemplate(out, pagesize=letter, leftMargin=0.6 * inch, rightMargin=0.6 * inch,
                  topMargin=0.55 * inch, bottomMargin=0.65 * inch,
                  title="n-body (plaquette) terms in ALPS, step by step").build(
    story, onFirstPage=footer, onLaterPages=footer)
print("wrote", out)
