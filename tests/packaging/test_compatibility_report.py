"""An expected difference must not hide a different or newly missing result."""

import importlib.util
from pathlib import Path

import pytest

path = Path(__file__).resolve().parents[2] / ".github/scripts/pyalps_compatibility/compare.py"
spec = importlib.util.spec_from_file_location("compatibility_compare", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def report(records, numpy="1.26.4"):
    return {"environment": {"python": "3.10.20", "numpy": numpy}, "records": records}


@pytest.mark.parametrize("candidate", [{"fixed": 3, "stable": 7}, {"fixed": 2}])
def test_classification_does_not_hide_changed_or_missing_values(candidate):
    expected = {"fixed": {"legacy": 1, "nanobind": 2, "reason": "known fix"}}
    result = module.compare(
        report({"fixed": 1, "stable": 7}), report(candidate), expected
    )
    assert result["unexpected"]


def test_comparison_accepts_only_the_classified_values():
    expected = {"fixed": {"legacy": 1, "nanobind": 2, "reason": "known fix"}}
    result = module.compare(
        report({"fixed": 1, "stable": 7}), report({"fixed": 2, "stable": 7}), expected
    )
    assert result["identical"] == 1 and result["classified_differences"] == 1
    assert not result["unexpected"] and not result["stale_classifications"]
    result = module.compare(report({"fixed": 2}), report({"fixed": 2}), expected)
    assert result["stale_classifications"] == ["fixed"]


def test_comparison_rejects_different_dependencies():
    with pytest.raises(ValueError, match="matching numpy"):
        module.compare(report({}), report({}, numpy="2.0.0"), {})
