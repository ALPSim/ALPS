import copy
import importlib
import sys
import numpy as np

legacy, action, filename = sys.argv[1] == "legacy", sys.argv[2], sys.argv[3]


def module(n):
    return importlib.import_module(n if legacy else "pyalps.cxx." + n)


module("pyalea_c")
d = module("pymcdata_c")
h = module("pyngshdf5_c")
p = module("pyngsparams_c")
r = module("pyngsrandom01_c")
o = module("pyngsobservables_c")
module("pyngsobservable_c")
results = module("pyngsresult_c")
ar = h.hdf5_archive_impl(filename, "w" if action == "write" else "r")
if action == "write":
    params = p.params(
        {
            "L": 16,
            "T": 1.25,
            "SEED": 42,
            "label": "check Ω",
            "couplings": np.array([1.0, 2.0, 3.0]),
        }
    )
    rng = r.random01(42)
    for _ in range(7):
        rng()
    expected_rng = copy.deepcopy(rng)
    expected = [expected_rng() for _ in range(8)]
    observables = o.observables()
    observables.createRealObservable("Energy")
    observables.createRealVectorObservable("Correlations")
    for i in range(64):
        observables["Energy"].append(float(i))
        observables["Correlations"].append(np.array([float(i), 2.0 * i]))
    result = results.observable2result(observables["Energy"])
    ar["/expected_rng"] = np.array(expected)
    for path, value in [
        ("parameters", params),
        ("rng", rng),
        ("observables", observables),
        ("result", result),
        ("vector_result", results.observable2result(observables["Correlations"])),
    ]:
        ar.set_context("/" + path)
        value.save(ar)
else:
    params = p.params(ar, "/parameters")
    assert params["L"] == 16 and params["T"] == 1.25 and params["label"] == "check Ω"
    np.testing.assert_array_equal(ar["/parameters/couplings"], [1, 2, 3])
    rng = r.random01()
    ar.set_context("/rng")
    rng.load(ar)
    np.testing.assert_array_equal([rng() for _ in range(8)], ar["/expected_rng"])
    observables = o.observables()
    observables.load(ar, "/observables")
    for name, mean in [("Energy", 31.5), ("Correlations", [31.5, 63.0])]:
        result = results.observable2result(observables[name])
        np.testing.assert_allclose(result.mean, mean)
        assert result.count == 64
    # Legacy result.load itself is broken. Inspect its persisted values, and
    # exercise actual scalar/vector restoration through the corrected binding.
    assert ar["/result/mean/value"] == 31.5 and ar["/result/count"] == 64
    if not legacy:
        for path, mean in [("/result", 31.5), ("/vector_result", [31.5, 63.0])]:
            restored = results.result()
            ar.set_context(path)
            restored.load(ar)
            np.testing.assert_allclose(restored.mean, mean)
            assert restored.count == 64
    print(("legacy" if legacy else "nanobind") + " read checkpoint: OK")
ar.close()
