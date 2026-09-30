import json
import os
import pathlib
import sys
import tempfile
import numpy as np
import pyalps

root = pathlib.Path(os.environ["PYALPS_WORKFLOW_ROOT"])
root.mkdir(parents=True, exist_ok=True)
app = sys.argv[1]
folder = pathlib.Path(tempfile.mkdtemp(prefix=f"{app}-", dir=root))
os.chdir(folder)
common = {
    "LATTICE": "chain lattice",
    "MODEL": "spin",
    "local_S": 0.5,
    "L": 4,
    "J": 1,
    "SEED": 42,
}
monte = {**common, "T": 1.0, "THERMALIZATION": 32, "SWEEPS": 256}
cases = {
    "loop": {**monte, "ALGORITHM": "loop"},
    # A short smoke run must not depend on thermalization growing the default
    # ten-slot operator string before measurement begins.
    "dirloop_sse": {**monte, "INITIAL_CUTOFF": 128},
    "spinmc": {
        "LATTICE": "square lattice",
        "MODEL": "Ising",
        "L": 4,
        "J": 1,
        "T": 2.0,
        "UPDATE": "cluster",
        "THERMALIZATION": 32,
        "SWEEPS": 256,
        "SEED": 42,
    },
    "sparsediag": {
        **common,
        "CONSERVED_QUANTUMNUMBERS": "Sz",
        "Sz_total": 0,
        "NUMBER_EIGENVALUES": 2,
    },
    "fulldiag": {
        **common,
        "CONSERVED_QUANTUMNUMBERS": "Sz",
        "T_MIN": 0.5,
        "T_MAX": 1.5,
        "DELTA_T": 0.5,
    },
    "dmrg": {
        **common,
        "LATTICE": "open chain lattice",
        "CONSERVED_QUANTUMNUMBERS": "N,Sz",
        "Sz_total": 0,
        "SWEEPS": 2,
        "MAXSTATES": 16,
        "NUMBER_EIGENVALUES": 1,
    },
}
input_file = pyalps.writeInputFiles("test", [cases[app]])
status = pyalps.runApplication(app, input_file, Tmin=1, T=45, writexml=True)
assert status[0] == 0, (app, status)
files = pyalps.getResultFiles(prefix="test")
assert files
if app in ("loop", "spinmc", "dirloop_sse"):
    data = pyalps.loadMeasurements(files, ["Energy"])
elif app in ("sparsediag", "fulldiag"):
    data = pyalps.loadSpectra(files)
else:
    data = pyalps.loadEigenstateMeasurements(files)
values = []
for dataset in pyalps.flatten(data):
    for item in np.asarray(dataset.y, dtype=object).flat:
        if isinstance(item, np.ndarray):
            values.extend(np.ravel(item).tolist())
        elif hasattr(item, "mean") and not callable(item.mean):
            values.append(item.mean)
        else:
            values.append(item)
assert values, (app, "no loaded data")
assert np.isfinite(np.asarray(values, dtype=float)).all(), (app, values)
if app in ("sparsediag", "fulldiag"):
    np.testing.assert_allclose(min(values), -2.0, rtol=0, atol=1e-10)
print(
    json.dumps(
        {
            "application": app,
            "files": len(files),
            "loaded_finite_values": len(values),
            "sample": list(map(float, values[:4])),
        }
    ),
    flush=True,
)
