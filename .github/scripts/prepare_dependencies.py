"""Prepare pinned CI dependencies using official sources and reusable caches.

CI builds selected Boost releases when the OS package is too old.
No project-owned GitHub release assets are required.
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


def prepare_boost(version: str, destination: Path):
    checksums = json.loads((ROOT / '.github/dependencies.json').read_text())['boost']
    if version not in checksums:
        raise ValueError(f'Unpinned Boost version: {version}')
    stamp = destination / '.alps-boost-sha256'
    fingerprint = checksums[version] + os.environ.get('CXX', 'c++') + os.environ.get('ALPS_BOOST_MPI', 'ON')
    if stamp.is_file() and stamp.read_text() == fingerprint:
        return
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(f'Existing dependency directory has a different configuration: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.alps-boost-', dir=destination.parent) as temporary:
        archive = Path(temporary) / 'boost.tar.gz'
        url = f'https://archives.boost.io/release/{version}/source/boost_{version.replace(".", "_")}.tar.gz'
        digest = hashlib.sha256()
        with urllib.request.urlopen(url, timeout=120) as response, archive.open('wb') as stream:
            while chunk := response.read(1024 * 1024):
                digest.update(chunk)
                stream.write(chunk)
        if digest.hexdigest() != checksums[version]:
            raise ValueError(f'Checksum mismatch for Boost {version}')
        source = Path(temporary) / 'source'
        with tarfile.open(archive) as package:
            package.extractall(source, filter='data')
        extracted, = source.iterdir()
        subprocess.run(['bash', str(ROOT / '.github/scripts/build-boost.sh'), str(extracted), str(destination.resolve())], check=True)
    stamp.write_text(fingerprint)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('package')
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    match = re.fullmatch(r'boost-(\d+\.\d+\.\d+)(?:-[a-z0-9-]+)?', args.package)
    if not match:
        parser.error('Expected boost-VERSION')
    prepare_boost(match[1], args.destination.resolve())


if __name__ == '__main__':
    main()
