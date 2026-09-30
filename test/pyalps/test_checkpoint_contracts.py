"""User workflows discovered by auditing the migration beyond API presence."""

import copy
import gc
import os
import subprocess
import sys
import textwrap
import weakref

import numpy as np
import pytest

from pyalps import hdf5, ngs


def result(vector=False, offset=0):
    observable = (ngs.createRealVectorObservable if vector else ngs.createRealObservable)("samples")
    for i in range(64):
        observable << (np.array([i + offset, 2.0 * i]) if vector else float(i + offset))
    return ngs.observable2result(observable)


def test_parameters_preserve_arithmetic_mutation_and_ownership(tmp_path):
    array = np.array([1.0, 2.0])
    reference = weakref.ref(array)
    parameters = ngs.params({"array": array, "list": [1, 2], "tuple": (3, 4)})
    np.testing.assert_array_equal(parameters["array"] * 2, [2, 4])
    array[0] = 9
    del array
    gc.collect()
    assert reference() is not None
    parameters["array"] += 1
    parameters["list"].append(3)
    parameters["list"][0] = 8
    assert parameters["tuple"] * 2 == (3, 4, 3, 4)
    duplicate = copy.deepcopy(parameters)
    duplicate["array"][0] = -1
    duplicate["list"].clear()
    np.testing.assert_array_equal(parameters["array"], [10, 3])
    assert parameters["list"] == [8, 2, 3]
    with hdf5.archive(str(tmp_path / "parameters.h5"), "w") as archive:
        archive["parameters"] = parameters
        np.testing.assert_array_equal(archive["parameters/array"], [10, 3])
        np.testing.assert_array_equal(archive["parameters/list"], [8, 2, 3])
    del parameters["array"]
    gc.collect()
    assert reference() is None


def test_parameter_checkpoint_retains_arrays_metadata_and_large_integers(tmp_path):
    values = {"array": np.arange(6.0).reshape(2, 3), "large": 2 ** 53 + 1,
              "metadata": {"label": "run", "ids": [1, 2]}, "flags": [True, False]}
    filename = str(tmp_path / "parameters.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["custom"] = ngs.params(values)
    with hdf5.archive(filename, "r") as archive:
        parameters = ngs.params(archive, "/custom")
        assert archive.context == "/"
    np.testing.assert_array_equal(parameters["array"] * 2, values["array"] * 2)
    assert parameters["large"] == 2 ** 53 + 1
    assert parameters["metadata"]["label"] == "run"
    assert parameters["flags"] == [True, False]
    parameters["array"][0, 0] = 99
    with hdf5.archive(filename, "a") as archive:
        archive["again"] = parameters
        assert archive["again/array"][0, 0] == 99


@pytest.mark.parametrize("vector", [False, True])
def test_result_loads_default_and_existing_objects_without_writing(tmp_path, vector):
    original = result(vector)
    filename = str(tmp_path / "result.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["result"] = original
    restored = ngs.result()
    assert restored.count == 0 and repr(restored) == "No Measurements"
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/result")
        restored.load(archive)
        assert archive.context == "/result"
    np.testing.assert_array_equal(restored.mean, original.mean)
    np.testing.assert_array_equal(restored.error, original.error)
    assert restored.count == original.count
    alias = ngs.result(restored)
    replacement = result(not vector, offset=100)
    with hdf5.archive(filename, "w") as archive:
        archive["result"] = replacement
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/result")
        restored.load(archive)
    np.testing.assert_array_equal(alias.mean, original.mean)
    np.testing.assert_array_equal(restored.mean, replacement.mean)


def test_result_collection_load_uses_context_and_keeps_old_references(tmp_path):
    first = ngs.results()
    first["energy/with&name"] = result()
    first["vector"] = result(True)
    filename = str(tmp_path / "results.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["custom/results"] = first
    restored = ngs.results()
    restored["old"] = result(offset=10)
    old_reference = restored["old"]
    with hdf5.archive(filename, "r") as archive:
        restored.load(archive, "/custom/results")
        assert archive.context == "/"
        assert set(restored) == set(first)
        energy_reference = restored["energy/with&name"]
        restored.load(archive, "/custom/results")
    del restored
    gc.collect()
    assert old_reference.mean == 41.5
    assert energy_reference.mean == 31.5


def test_failed_loads_preserve_values_and_context(tmp_path):
    filename = str(tmp_path / "broken.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["bad/mean/value"] = 100.0  # no mandatory count dataset
    value = result()
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/bad")
        with pytest.raises(Exception):
            value.load(archive)
        assert archive.context == "/bad"
        assert value.mean == 31.5
        archive.set_context("/")
        parameters = ngs.params({"kept": [1, 2]})
        with pytest.raises(Exception):
            parameters.load(archive, "/missing")
        assert archive.context == "/" and parameters["kept"] == [1, 2]
        results = ngs.results()
        results["kept"] = value
        with pytest.raises(Exception):
            results.load(archive, "/bad")
        assert archive.context == "/" and results["kept"].mean == 31.5


@pytest.mark.parametrize("vector", [False, True])
def test_result_unary_operations_and_deepcopy_do_not_mutate_original(vector):
    original = result(vector)
    mean = np.array(original.mean)
    positive, negative, duplicate = +original, -original, copy.deepcopy(original)
    duplicate += 10
    np.testing.assert_array_equal(original.mean, mean)
    np.testing.assert_array_equal(positive.mean, mean)
    np.testing.assert_array_equal(negative.mean, -mean)
    np.testing.assert_array_equal(duplicate.mean, mean + 10)


def test_empty_result_operations_raise_instead_of_crashing():
    code = '''
        from pyalps import ngs
        result = ngs.result()
        for operation in (lambda: result + 1, lambda: abs(result),
                          lambda: result ** 2, lambda: result.sin(),
                          lambda: result.mean, lambda: result.error):
            try:
                operation()
            except (RuntimeError, ValueError):
                pass
            else:
                raise AssertionError("empty result operation should fail")
    '''
    completed = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", textwrap.dedent(code)],
        capture_output=True, text=True, timeout=30,
        env={**os.environ, "MallocScribble": "1"},
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
