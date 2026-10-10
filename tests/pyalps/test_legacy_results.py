"""Read frozen application output independently of the current SDK writer."""
import gzip
import json
from pathlib import Path

import numpy as np
import pytest
import pyalps

FIXTURES = Path(__file__).with_name("fixtures") / "application_results"
REFERENCES = json.loads((FIXTURES / "expected.json").read_text())


@pytest.mark.parametrize("application", REFERENCES)
def test_legacy_application_results(application, tmp_path):
    archive = tmp_path / f"{application}.h5"
    archive.write_bytes(gzip.decompress((FIXTURES / f"{application}.h5.gz").read_bytes()))
    monte_carlo = application in ("spinmc", "looper")
    loader = pyalps.loadMeasurements if monte_carlo else pyalps.loadEigenstateMeasurements
    data = pyalps.flatten(loader([str(archive)]))
    assert data, "The loader silently discarded historical results"
    assert all(dataset.props["L"] == 4 for dataset in data)
    energy = [dataset for dataset in data if dataset.props["observable"] == "Energy"]
    assert energy
    expected = REFERENCES[application]
    if monte_carlo:
        assert len(energy) == 1
        value = energy[0].y[0]
        np.testing.assert_allclose([value.mean, value.error],
                                   [expected["mean"], expected["error"]], rtol=1e-14)
        assert value.count == expected["count"]
        np.testing.assert_array_equal(energy[0].x, [0])
    else:
        values = sorted(float(value) for dataset in energy for value in dataset.y)
        np.testing.assert_allclose(values, expected["energies"], rtol=1e-14, atol=1e-15)
        if application == "dmrg":
            local = next(dataset for dataset in data
                         if dataset.props["observable"] == "Local magnetization")
            np.testing.assert_array_equal(local.x, np.arange(4))
            np.testing.assert_allclose(local.y, expected["local"], rtol=1e-14, atol=1e-15)
