"""Check source integrity and cache reuse without downloading packages."""
import importlib.util
import io
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('prepare_dependencies', ROOT / '.github/scripts/prepare_dependencies.py')
dependencies = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dependencies)


def test_boost_checksum_failure_never_extracts_or_builds(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'corrupt archive'))
    def forbidden(*a, **kw):
        pytest.fail('An unverified dependency must not be extracted or built')
    monkeypatch.setattr(dependencies.tarfile, 'open', forbidden)
    monkeypatch.setattr(dependencies.subprocess, 'run', forbidden)
    with pytest.raises(ValueError, match='Checksum mismatch'):
        dependencies.prepare_boost('1.91.0', tmp_path / 'installed')
    assert not (tmp_path / 'installed').exists()


def test_matching_boost_cache_avoids_network(tmp_path, monkeypatch):
    monkeypatch.setenv('CXX', 'clang++')
    monkeypatch.setenv('ALPS_BOOST_MPI', 'OFF')
    checksum = json.loads((ROOT / '.github/dependencies.json').read_text())['boost']['1.91.0']
    (tmp_path / '.alps-boost-sha256').write_text(checksum + 'clang++OFF')
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached dependency downloaded again'))
    dependencies.prepare_boost('1.91.0', tmp_path)
