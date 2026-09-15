"""Check release versions before building or publishing (Python >= 3.11)."""

import argparse
from email.parser import BytesParser
import os
from pathlib import Path
import re
import tarfile
import tomllib
import zipfile

from packaging.utils import (
    canonicalize_name,
    parse_sdist_filename,
    parse_wheel_filename,
)
from packaging.version import Version


CORE_PATTERN = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
TAG_PATTERN = rf"v{CORE_PATTERN}(?:-(?:alpha|beta|rc|dev)\.[0-9]+)?"


def check_version(root: Path, ref: str) -> Version:
    core = (root / "ALPS_VERSION.txt").read_text().strip()
    if not re.fullmatch(CORE_PATTERN, core):
        raise ValueError("ALPS_VERSION.txt must contain MAJOR.MINOR.PATCH")

    with (root / "pyproject.toml").open("rb") as stream:
        project = tomllib.load(stream)["project"]
    version = Version(project["version"])
    if version.release != Version(core).release:
        raise ValueError(
            f"pyproject.toml version {version} disagrees with ALPS_VERSION.txt ({core})"
        )

    if ref.startswith("refs/tags/"):
        tag = ref.removeprefix("refs/tags/")
        if not re.fullmatch(TAG_PATTERN, tag):
            raise ValueError(
                f"Invalid release tag {tag!r}; expected vMAJOR.MINOR.PATCH"
                "[-{alpha,beta,rc,dev}.N]"
            )
        if Version(tag) != version:
            raise ValueError(
                f"Release tag {tag} disagrees with pyproject.toml ({version})"
            )

    return version


def check_distributions(directory: Path, expected: Version) -> None:
    artifacts = sorted(directory.iterdir())
    if not artifacts:
        raise ValueError(f"No distributions found in {directory}")

    for artifact in artifacts:
        if artifact.suffix == ".whl":
            name, version, _, _ = parse_wheel_filename(artifact.name)
            with zipfile.ZipFile(artifact) as archive:
                members = [
                    member for member in archive.namelist()
                    if member.count("/") == 1 and member.endswith(".dist-info/METADATA")
                ]
                if len(members) != 1:
                    raise ValueError(f"{artifact.name}: expected one wheel METADATA file")
                metadata = archive.read(members[0])
        elif artifact.name.endswith(".tar.gz"):
            name, version = parse_sdist_filename(artifact.name)
            with tarfile.open(artifact) as archive:
                members = [
                    member for member in archive.getmembers()
                    if member.isfile() and member.name.count("/") == 1
                    and member.name.endswith("/PKG-INFO")
                ]
                if len(members) != 1:
                    raise ValueError(f"{artifact.name}: expected one sdist PKG-INFO file")
                with archive.extractfile(members[0]) as stream:
                    metadata = stream.read()
        else:
            raise ValueError(f"Unexpected distribution: {artifact.name}")

        if name != "pyalps" or version != expected:
            raise ValueError(f"{artifact.name}: expected a pyalps {expected} distribution")
        headers = BytesParser().parsebytes(metadata)
        if (
            canonicalize_name(headers.get("Name", "")) != "pyalps"
            or Version(headers.get("Version", "")) != expected
        ):
            raise ValueError(f"{artifact.name}: metadata does not describe pyalps {expected}")
        print(f"Verified {artifact.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--ref", default=os.environ.get("GITHUB_REF", ""))
    parser.add_argument("--dist", type=Path)
    args = parser.parse_args()
    try:
        version = check_version(args.root, args.ref)
        print(f"Release version: {version}")
        if args.dist is not None:
            check_distributions(args.dist, version)
    except (OSError, ValueError, KeyError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Release version check failed: {error}\n")


if __name__ == "__main__":
    main()
