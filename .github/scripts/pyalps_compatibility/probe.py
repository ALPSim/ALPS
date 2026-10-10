"""Run the same public binding probes against separately built old/new modules."""

import copy
import importlib
import hashlib
import platform
import json
import os
import sys
import tempfile
import warnings
import numpy as np

# Execute old/new in separate processes with matching Python and NumPy.

legacy = sys.argv[1] == "legacy"


def module(name):
    return importlib.import_module(name if legacy else "pyalps.cxx." + name)


a = module("pyalea_c")
d = module("pymcdata_c")
h = module("pyngshdf5_c")
p = module("pyngsparams_c")
t = module("pytools_c")
r = module("pyngsrandom01_c")
results = {}


def normalize(v):
    if isinstance(v, np.ndarray):
        return {
            "array": normalize(v.tolist()),
            "dtype": v.dtype.str,
            "shape": list(v.shape),
        }
    if isinstance(v, np.generic):
        return normalize(v.item())
    if isinstance(v, float):
        if not np.isfinite(v):
            return str(v)
        return float(f"{v:.11g}")
    if isinstance(v, complex):
        return {"real": normalize(v.real), "imag": normalize(v.imag)}
    if isinstance(v, (str, int, bool)) or v is None:
        return v
    if isinstance(v, dict):
        return {str(k): normalize(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [normalize(x) for x in v]
    if hasattr(v, "timeseries"):
        return normalize(v.timeseries())
    if hasattr(v, "first"):
        return [normalize(v.first), normalize(v.second)]
    if hasattr(v, "mean"):
        return {"mean": normalize(v.mean), "error": normalize(v.error)}
    return str(v)


def probe(name, call):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results[name] = {"ok": normalize(call())}
    except Exception as e:
        results[name] = {"error": type(e).__name__, "message": str(e).splitlines()[0]}


for name in (
    "pyalea_c",
    "pymcdata_c",
    "pytools_c",
    "pyngsparams_c",
    "pyngshdf5_c",
    "pyngsbase_c",
    "pyngsobservable_c",
    "pyngsobservables_c",
    "pyngsresult_c",
    "pyngsresults_c",
    "pyngsapi_c",
    "pyngsrandom01_c",
    "pyngsaccumulator_c",
):
    m = module(name)
    results["surface/" + name] = {
        "ok": sorted(n for n in dir(m) if not n.startswith("_"))
    }
    for n in dir(m):
        cl = getattr(m, n)
        if isinstance(cl, type) and not n.startswith("_"):
            results["members/" + name + "/" + n] = {
                "ok": sorted(x for x in dir(cl) if not x.startswith("_"))
            }
for cls, mean, error in [
    (d.MCScalarData, 1.5, 0.25),
    (d.MCVectorData, np.array([1.5, 2.5, 3.5]), np.array([0.25, 0.5, 0.75])),
]:
    value = cls(mean, error)
    name = cls.__name__
    for method in (
        "sq",
        "cb",
        "sqrt",
        "cbrt",
        "exp",
        "log",
        "sin",
        "cos",
        "tan",
        "sinh",
        "cosh",
        "tanh",
    ):
        probe(name + "/" + method, lambda method=method: getattr(value, method)())
    for othername, other in [
        ("self", value),
        ("float", 2.0),
        ("int", 2),
        ("npfloat", np.float64(2.0)),
    ]:
        for method in (
            "__add__",
            "__sub__",
            "__mul__",
            "__truediv__",
            "__radd__",
            "__rsub__",
            "__rmul__",
            "__rtruediv__",
        ):
            probe(
                name + "/" + method + "/" + othername,
                lambda method=method, other=other: getattr(value, method)(other),
            )
    for label, fn in [
        ("neg", lambda: -value),
        ("pos", lambda: +value),
        ("abs", lambda: abs(value)),
        ("pow", lambda: value**2),
        ("copy", lambda: copy.deepcopy(value)),
        ("repr", lambda: repr(value)),
        ("format", lambda: format(value, ".3f")),
    ]:
        probe(name + "/" + label, fn)
    if cls is d.MCVectorData:
        for label, index in [
            ("first", 0),
            ("last", -1),
            ("slice", slice(0, 2)),
            ("full", slice(None)),
            ("empty", slice(2, 1)),
        ]:
            probe(name + "/" + label, lambda index=index: value[index])
for cls, vector in [
    (a.RealObservable, False),
    (a.RealTimeSeriesObservable, False),
    (a.RealVectorObservable, True),
    (a.RealVectorTimeSeriesObservable, True),
]:
    obs = cls("samples")
    for i in range(256):
        x = float(np.sin(i * 0.73) + 2)
        obs << (np.array([x, 2 * x]) if vector else x)
    for key in ("mean", "error", "tau", "variance", "count", "converged_errors"):
        probe(cls.__name__ + "/" + key, lambda key=key: getattr(obs, key))
    with tempfile.TemporaryDirectory() as directory:
        filename = os.path.join(directory, "sample.h5")
        obs.save(filename)
        data = d.MCVectorData() if vector else d.MCScalarData()
        data.load(filename, "/simulation/results/samples")
        cls_ts = a.MCVectorTimeseries if vector else a.MCScalarTimeseries
        probe(cls.__name__ + "/to_timeseries", lambda: cls_ts(data))
series = a.MCScalarTimeseries(np.sin(np.arange(256) * 0.73) + 2)
for f in (
    "mean",
    "variance",
    "uncorrelated_error",
    "binning_error",
    "running_mean",
    "reverse_running_mean",
):
    probe("timeseries/" + f, lambda f=f: getattr(a, f)(series))
for f, arg in [
    ("autocorrelation_distance", 20),
    ("cut_head_distance", 10),
    ("cut_tail_distance", 10),
    ("cut_head_limit", 0.5),
    ("cut_tail_limit", 0.5),
]:
    probe("timeseries/" + f, lambda f=f, arg=arg: getattr(a, f)(series, arg))
for cls in (t.rng, r.random01):
    for seed in (0, 1, 42, 2147483647):
        rng = cls(seed)
        probe(
            "rng/" + cls.__name__ + "/" + str(seed), lambda: [rng() for _ in range(25)]
        )
for value in ("normal", "with / slash", "unicode Ω/漢字", "percent%", "@attribute", ""):
    probe(
        "name/" + value,
        lambda value=value: t.hdf5_name_decode(t.hdf5_name_encode(value)),
    )

with tempfile.TemporaryDirectory() as directory:
    archive = h.hdf5_archive_impl(os.path.join(directory, "data.h5"), "w")

    def roundtrip(value):
        archive["/value"] = value
        return archive["/value"]

    scalar_values = [
        True,
        False,
        0,
        -13,
        4.5,
        2 + 3j,
        "unicode Ω/漢字",
        "",
        np.int8(-7),
        np.uint8(240),
        np.int16(-999),
        np.uint16(60000),
        np.int32(-123456),
        np.uint32(3000000000),
        np.int64(-(2**40)),
        np.uint64(2**63 + 1),
        np.float32(1.25),
        np.float64(1.25),
        np.complex64(2 + 3j),
        np.complex128(2 + 3j),
    ]
    for i, value in enumerate(scalar_values):
        probe("hdf/scalar/" + str(i), lambda value=value: roundtrip(value))
    for dtype in (
        "bool",
        "int8",
        "uint8",
        "int16",
        "uint16",
        "int32",
        "uint32",
        "int64",
        "uint64",
        "float32",
        "float64",
        "complex64",
        "complex128",
    ):
        for shape in ((), (6,), (2, 3), (2, 1, 3)):
            data = np.arange(int(np.prod(shape))).astype(dtype).reshape(shape)
            for layout in ("C", "F", "readonly"):
                value = (
                    np.asfortranarray(data)
                    if layout == "F" and len(shape) > 1
                    else data.copy()
                )
                if layout == "readonly":
                    value.flags.writeable = False
                probe(
                    "hdf/array/" + dtype + "/" + str(shape) + "/" + layout,
                    lambda value=value: roundtrip(value),
                )
    for label, value in [
        ("list", [1, 2, 3]),
        ("nested", [[1.0, 2.0], [3.0, 4.0]]),
        ("ragged", [[1], [2, 3]]),
        ("mixed", [1, "two", 3.0]),
        ("dictionary", {"a": [1, 2], "b": {"c": 1.5}}),
        ("empty-list", []),
        ("empty-dict", {}),
        ("complex-list", [1 + 2j, 3 + 4j]),
    ]:
        probe("hdf/container/" + label, lambda value=value: roundtrip(value))
    for label, value in [
        ("python-numpy-integer-rows", [[1, 2], [np.int32(3), np.int32(4)]]),
        ("python-numpy-float-rows", [[1.0, 2.0], [np.float64(3.0), np.float64(4.0)]]),
        ("list-ndarray-dtype", [[1, 2], np.array([3, 4], dtype=np.int32)]),
    ]:
        probe("hdf/container/" + label, lambda value=value: roundtrip(value))
    archive.close()

for label, value in [
    ("int", 7),
    ("large-int", 2**40),
    ("float", 1.25),
    ("bool", True),
    ("str", "a"),
    ("complex", 1 + 2j),
    ("list", [1.0, 2.0]),
    ("array", np.array([1.0, 2.0])),
    ("matrix", np.ones((2, 2))),
    ("None", None),
    ("dict", {"x": 1}),
    ("tuple", (1, 2)),
    ("numpy-scalar", np.int64(7)),
]:

    def params_probe(value=value):
        par = p.params({"x": value})
        got = par["x"]
        return {"type": type(got).__name__, "value": normalize(got)}

    probe("params/" + label, params_probe)
modules = {}
for name in ("pyalea_c", "pymcdata_c", "pyngshdf5_c", "pyngsparams_c"):
    filename = module(name).__file__
    with open(filename, "rb") as source:
        modules[name] = {
            "path": filename,
            "sha256": hashlib.sha256(source.read()).hexdigest(),
        }
environment = {
    "python": platform.python_version(),
    "numpy": np.__version__,
    "platform": platform.platform(),
    "modules": modules,
}
with open(sys.argv[2], "w") as output:
    json.dump(
        {"environment": environment, "records": results},
        output,
        indent=2,
        sort_keys=True,
    )
print(
    len(results),
    "probes",
    len([v for v in results.values() if "error" in v]),
    "exceptions",
)
