"""Prepare pinned CI dependencies using official sources and reusable caches.

CI installs the headers of selected Boost releases when the OS package is too old,
and builds NFFT3 (for CT-HYB-SEGMENT's ALPS_ENABLE_NFFT), which most distributions
do not package. No project-owned GitHub release assets are required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]


def checksums(package: str) -> dict:
    return json.loads((ROOT / '.github/dependencies.json').read_text())[package]


def download_and_extract(url: str, checksum: str, temporary: Path, name: str) -> Path:
    """Download url, verify its SHA-256 before anything reads it, and extract it."""
    archive = temporary / 'archive.tar.gz'
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as response, archive.open('wb') as stream:
        while chunk := response.read(1024 * 1024):
            digest.update(chunk)
            stream.write(chunk)
    if digest.hexdigest() != checksum:
        raise ValueError(f'Checksum mismatch for {name}')
    source = temporary / 'source'
    with tarfile.open(archive) as package:
        if hasattr(tarfile, 'data_filter'):
            package.extractall(source, filter='data')
        else:
            # Python without extraction filters (Debian 12's 3.11.2): the archive's
            # checksum was verified against the pinned value above.
            package.extractall(source)  # noqa: S202
    extracted, = source.iterdir()
    return extracted


def prepare(package: str, version: str, destination: Path, url: str, build):
    """Build a pinned package into destination once; a stamp makes reruns free."""
    pinned = checksums(package)
    if version not in pinned:
        raise ValueError(f'Unpinned {package} version: {version}')
    stamp = destination / f'.alps-{package}-sha256'
    fingerprint = pinned[version] + os.environ.get('CXX', 'c++')
    if stamp.is_file() and stamp.read_text() == fingerprint:
        return
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(f'Existing dependency directory has a different configuration: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.alps-{package}-', dir=destination.parent) as temporary:
        source = download_and_extract(url, pinned[version], Path(temporary), f'{package} {version}')
        build(source, destination.resolve())
    stamp.write_text(fingerprint)


def prepare_boost(version: str, destination: Path):
    def build(source, prefix):
        # ALPS uses Boost's headers and CMake package files, and no compiled Boost library.
        subprocess.run(['./bootstrap.sh'], cwd=source, check=True)
        subprocess.run(['./b2', '-d0', f'--prefix={prefix}', '--with-headers', 'install'], cwd=source, check=True)
    url = f'https://archives.boost.io/release/{version}/source/boost_{version.replace(".", "_")}.tar.gz'
    prepare('boost', version, destination, url, build)


def prepare_nfft(version: str, destination: Path):
    def build(source, prefix):
        # The release archive carries its configure script. Only the library is
        # needed; FFTW3 comes from the distribution.
        subprocess.run(['./configure', f'--prefix={prefix}', '--disable-examples',
                        '--disable-applications', '--disable-julia'], cwd=source, check=True)
        subprocess.run(['make', f'-j{os.cpu_count() or 1}'], cwd=source, check=True)
        subprocess.run(['make', 'install'], cwd=source, check=True)
    url = f'https://github.com/NFFT/nfft/releases/download/{version}/nfft-{version}.tar.gz'
    prepare('nfft', version, destination, url, build)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('package')
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    match = re.fullmatch(r'(boost|nfft)-(\d+\.\d+\.\d+)(?:-[a-z0-9-]+)?', args.package)
    if not match:
        parser.error('Expected boost-VERSION or nfft-VERSION')
    {'boost': prepare_boost, 'nfft': prepare_nfft}[match[1]](match[2], args.destination.resolve())


if __name__ == '__main__':
    main()
