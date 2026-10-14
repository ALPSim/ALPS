"""Run ALPS ED on a lattice + custom Hamiltonian, return sorted energies."""

import os
import re
import shutil
import subprocess
import tempfile

# repository root: this file is src/alps/model/docs/n_body_terms/benchmarks/alps_checks/alpsrun.py
ALPS = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 7))
BIN = os.environ.get("ALPS_BIN", f"{ALPS}/_build/default/bin")
MODELS = f"{ALPS}/src/alps/resources/models.xml"

os.environ.setdefault("ALPS_XML_PATH", f"{ALPS}/src/alps/resources")

_E = re.compile(r'<SCALAR_AVERAGE name="Energy"><MEAN>([^<]+)</MEAN>')


def ladder_lattice(name="ms-ladder", plaquette=True, boundary="open"):
    plaq = """
    <PLAQUETTE type="0">
      <SITE vertex="1" offset="0"/>
      <SITE vertex="2" offset="0"/>
      <SITE vertex="2" offset="1"/>
      <SITE vertex="1" offset="1"/>
    </PLAQUETTE>""" if plaquette else ""
    return f"""<LATTICES>
  <LATTICE name="chain lattice" dimension="1">
    <PARAMETER name="a" default="1"/>
    <BASIS><VECTOR>a</VECTOR></BASIS>
    <RECIPROCALBASIS><VECTOR>2*pi/a</VECTOR></RECIPROCALBASIS>
  </LATTICE>
  <UNITCELL name="ms-ladder-cell" dimension="1" vertices="2">
    <VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX>
    <VERTEX id="2"><COORDINATE>0.5</COORDINATE></VERTEX>
    <EDGE type="1"><SOURCE vertex="1" offset="0"/><TARGET vertex="2" offset="0"/></EDGE>
    <EDGE type="0"><SOURCE vertex="1" offset="0"/><TARGET vertex="1" offset="1"/></EDGE>
    <EDGE type="0"><SOURCE vertex="2" offset="0"/><TARGET vertex="2" offset="1"/></EDGE>{plaq}
  </UNITCELL>
  <LATTICEGRAPH name="{name}">
    <FINITELATTICE>
      <LATTICE ref="chain lattice"/>
      <EXTENT dimension="1" size="L"/>
      <BOUNDARY type="{boundary}"/>
    </FINITELATTICE>
    <UNITCELL ref="ms-ladder-cell"/>
  </LATTICEGRAPH>
</LATTICES>
"""


def model_library(hamiltonian_xml):
    body = open(MODELS).read()
    body = body[body.index("<MODELS>") + len("<MODELS>"):body.rindex("</MODELS>")]
    return "<MODELS>\n" + body + "\n" + hamiltonian_xml + "\n</MODELS>\n"


def run(lattice_xml, hamiltonian_xml, params, program="fulldiag", keep=False):
    tmp = tempfile.mkdtemp(prefix="ms_")
    try:
        with open(os.path.join(tmp, "lattices.xml"), "w") as f:
            f.write(lattice_xml)
        if hamiltonian_xml is not None:
            with open(os.path.join(tmp, "models.xml"), "w") as f:
                f.write(model_library(hamiltonian_xml))
        p = {"LATTICE_LIBRARY": "lattices.xml"}
        if hamiltonian_xml is not None:
            p["MODEL_LIBRARY"] = "models.xml"
        p.update(params)
        with open(os.path.join(tmp, "parm"), "w") as f:
            for k, v in p.items():
                f.write(f'{k}="{v}"\n' if isinstance(v, str) else f"{k}={v!r}\n")
            f.write("{ }\n")
        r = lambda *a: subprocess.run(a, cwd=tmp, check=True, capture_output=True, text=True)
        r(os.path.join(BIN, "parameter2xml"), "parm")
        try:
            out = r(os.path.join(BIN, program), "--write-xml", "parm.in.xml")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(e.stderr + e.stdout) from None
        files = [f for f in os.listdir(tmp) if f.endswith(".out.xml") and ".task" in f]
        if not files:
            raise RuntimeError("no task output:\n" + out.stdout[-800:] + out.stderr[-800:])
        xml = open(os.path.join(tmp, files[0])).read()
        return sorted(float(x) for x in _E.findall(xml))
    finally:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)
