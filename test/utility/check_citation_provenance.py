"""Run native HDF5 regressions in a temporary directory; no persistent test data."""
import subprocess
import sys
import tempfile
import importlib.util
from pathlib import Path


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


with tempfile.TemporaryDirectory(prefix="alps-citation-provenance-") as directory:
    subprocess.run([sys.argv[1], directory], check=True, timeout=60)
    # Optional independent HDF5 reader; no additional runtime dependency in pyalps.
    try:
        import h5py
    except ImportError:
        h5py = None
    if h5py is not None:
        root = Path(__file__).resolve().parents[2]
        reader = load("saved_citations", root / "lib/pyalps/citations.py")
        generator = load("citation_generator", root / "script/generate_citations.py")
        cff, policy, references, framework = generator.load_catalog(root)
        for filename in Path(directory).glob("*.h5"):
            with h5py.File(filename, "r") as source:
                if filename.name in ("future.h5", "fractional-schema.h5") or filename.name.startswith("invalid-complete-"):
                    try:
                        reader.read_citations(source)
                    except ValueError:
                        continue
                    raise AssertionError("Invalid metadata accepted: " + filename.name)
                for record in reader.read_citations(source):
                    expected = generator.snapshot(cff, policy, references, framework,
                                                  record["component"], record["software_version"], record["activity"])
                    assert record == expected, filename
                    path = "/provenance/alps/citations/records/" + record["id"]
                    for field in ("bibliography_cff", "notice"):
                        assert h5py.check_string_dtype(source[path + "/" + field].dtype).encoding == "utf-8"
        with h5py.File(Path(directory) / "copied.h5", "r") as source:
            record = reader.read_citations(source)[0]
            path = Path(directory) / "CITATION.cff"
            reader.export_citations(source, path, record["id"])
            assert path.read_text(encoding="utf-8") == record["bibliography_cff"]
