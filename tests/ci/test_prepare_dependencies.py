"""Check source integrity and architecture selection without downloading packages."""
import hashlib
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
    checksum = json.loads((ROOT / '.github/ci-matrix.json').read_text())['boost']['1.91.0']
    (tmp_path / '.alps-boost-sha256').write_text(checksum + 'clang++OFF')
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached dependency downloaded again'))
    dependencies.prepare_boost('1.91.0', tmp_path)


@pytest.mark.parametrize('architecture', ['x64', 'arm64'])
def test_windows_uses_pinned_baseline_and_matching_triplet(tmp_path, monkeypatch, architecture):
    commands = []
    monkeypatch.setattr(dependencies.subprocess, 'run', lambda args, **kw: commands.append(args))
    dependencies.prepare_windows(architecture, tmp_path)
    baseline = json.loads((ROOT / 'cmake/vcpkg/vcpkg.json').read_text())['builtin-baseline']
    assert ['git', '-C', str(tmp_path), 'checkout', '--detach', baseline] in commands
    install = commands[-1]
    assert f'--triplet={architecture}-windows' in install
    assert f'--x-manifest-root={ROOT / "cmake/vcpkg"}' in install
    assert '--x-feature=tests' in install
    assert any(arg.startswith('--overlay-ports=') for arg in install) == (architecture == 'arm64')
