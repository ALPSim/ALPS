# SPDX-License-Identifier: MIT
"""Read historical ALPS citation snapshots without consulting today's catalog.

Paths use the existing ALPS HDF5 bindings. An already-open ALPS archive or h5py
file is also accepted. No YAML parser, network access, or write operation is used.
"""
import hashlib
import json
from numbers import Integral
import os
from pathlib import Path
import re

_ROOT = "/provenance/alps/citations"
_TEXT_FIELDS = ("component", "software_version", "catalog_sha256", "selection_basis",
                "activity", "bibliography_cff", "notice", "request", "note", "license_note")
_ROLES = ("algorithm", "implementation", "framework")


def _text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if not isinstance(value, str):
        raise ValueError("Citation text must be a UTF-8 string")
    return str(value)


def _value(archive, path):
    value = archive[path]
    # h5py returns datasets; ALPS returns the value directly.
    return value[()] if hasattr(value, "shape") and hasattr(value, "file") else value


def _group(archive, path):
    if hasattr(archive, "is_group"):
        return archive.is_group(path)
    return path in archive and hasattr(archive[path], "keys")


def _data(archive, path):
    if hasattr(archive, "is_data"):
        return archive.is_data(path)
    return path in archive and not _group(archive, path)


def _children(archive, path):
    return archive.list_children(path) if hasattr(archive, "list_children") else archive[path].keys()


def read_citations(source):
    """Return validated snapshot dictionaries, sorted by stable id.

    Old files without citation metadata return []. Unsupported schemas and
    damaged snapshots raise ValueError. Reading never updates the saved guidance.
    Accept a filename, an open pyalps.hdf5.archive, or an open h5py.File.
    """
    if isinstance(source, (str, bytes, os.PathLike)):
        from .hdf5 import archive
        with archive(os.fsdecode(source), "r") as handle:
            return read_citations(handle)
    try:
        return _read_citations(source)
    except (LookupError, TypeError, UnicodeError) as error:
        raise ValueError("Malformed ALPS citation provenance: " + str(error)) from error


def _read_citations(source):
    if not _group(source, _ROOT):
        if _data(source, _ROOT):
            raise ValueError("Citation provenance path is not a group")
        return []
    version = _value(source, _ROOT + "/schema_version")
    if not isinstance(version, Integral):
        raise ValueError("Citation schema version must be a scalar integer")
    if version != 1:
        raise ValueError(f"Unsupported ALPS citation provenance schema: {version}")
    records_path = _ROOT + "/records"
    if not _group(source, records_path):
        if _data(source, records_path):
            raise ValueError("Citation records path is not a group")
        return []
    records = []
    for identity in sorted(_children(source, records_path)):
        path = records_path + "/" + identity
        if not _group(source, path):
            raise ValueError("Citation snapshot path is not a group")
        if not _data(source, path + "/complete"):
            continue
        complete = _value(source, path + "/complete")
        if not isinstance(complete, Integral) or complete not in (0, 1):
            raise ValueError("Invalid citation completion marker")
        if complete == 0:
            continue
        record = {field: _text(_value(source, path + "/" + field)) for field in _TEXT_FIELDS}
        for role in _ROLES:
            role_path = path + "/" + role
            # ALPS represents empty vectors with the HDF5 null dataspace.
            if (hasattr(source, "is_null") and source.is_null(role_path)) or (
                    not hasattr(source, "is_null") and source[role_path].shape is None):
                values = []
            else:
                values = _value(source, role_path)
            record[role] = [_text(value) for value in values]
        payload = json.dumps(record, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        if not re.fullmatch(r"[0-9a-f]{64}", identity) or hashlib.sha256(payload).hexdigest() != identity:
            raise ValueError("Citation snapshot fingerprint mismatch: " + identity)
        if not re.fullmatch(r"[0-9a-f]{64}", record["catalog_sha256"]):
            raise ValueError("Invalid citation catalog fingerprint")
        if not all(record[field] for field in ("component", "software_version", "bibliography_cff", "notice", "framework")):
            raise ValueError("Citation snapshot has an empty required field")
        if record["selection_basis"] != "application" or record["activity"] not in ("calculation", "analysis", "unspecified"):
            raise ValueError("Invalid citation selection/activity")
        record["id"] = identity
        records.append(record)
    return records


def export_citations(source, destination, snapshot_id=None):
    """Export the saved bibliography as a standard CITATION.cff document.

    A file with different bibliographies requires an explicit snapshot_id;
    identical calculation/analysis bibliographies can share a single export.
    Returns the destination Path. Existing destinations are not overwritten.
    """
    records = read_citations(source)
    if snapshot_id is not None:
        records = [record for record in records if record["id"] == snapshot_id]
        if not records:
            raise ValueError("Unknown citation snapshot: " + snapshot_id)
    documents = {record["bibliography_cff"] for record in records}
    if not documents:
        raise ValueError("No saved citation information is available")
    if len(documents) != 1:
        raise ValueError("Several saved bibliographies are available; select a snapshot_id")
    path = Path(destination)
    with path.open("x", encoding="utf-8", newline="\n") as output:
        output.write(documents.pop())
    return path
