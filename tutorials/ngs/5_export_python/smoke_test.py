#!/usr/bin/env python3
"""Exercise the public downstream simulation-export compatibility helper."""

import os
import tempfile

# Importing the consumer first verifies its wheel-runtime rpath. Its module
# initializer loads the owning pyalps bindings before registering C++ types.
import ising_c
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs


parameters = ngs.params({"SEED": 7, "SWEEPS": 10})
simulation = ising_c.sim(parameters)

assert issubclass(ising_c.sim, ngs.mcbase)
assert isinstance(simulation, ngs.mcbase)
assert int(simulation.parameters["SWEEPS"]) == 10
assert len(simulation.measurements) == 1
assert 0.0 <= simulation.random() < 1.0
# The base descriptors must work for a downstream C++ simulation too. This
# used to throw std::bad_cast because mcbase assumed every instance was its
# Python trampoline alias.
assert ngs.mcbase.parameters.__get__(simulation) is simulation.parameters
assert ngs.mcbase.measurements.__get__(simulation) is simulation.measurements
assert ngs.mcbase.random.__get__(simulation) is simulation.random
assert simulation.run(lambda: False)
assert simulation.resultNames() == ["Magnetization"]
before = simulation.collectResults()
assert before["Magnetization"].count == 10

with tempfile.TemporaryDirectory() as directory:
    checkpoint = os.path.join(directory, "ising.h5")
    with hdf5.archive(checkpoint, "w") as archive:
        simulation.save(archive)

    restored = ising_c.sim(parameters)
    with hdf5.archive(checkpoint, "r") as archive:
        restored.load(archive)

    after = restored.collectResults()
    assert restored.resultNames() == simulation.resultNames()
    assert after["Magnetization"].count == before["Magnetization"].count
    assert after["Magnetization"].mean == before["Magnetization"].mean

    # `archive[path] = simulation` must write exactly what simulation.save()
    # writes. It reaches save() through the archive-savable marker that this
    # class inherits from mcbase -- but mcbase's own bound save() is a
    # deliberately *non-virtual* qualified call (so that a Python subclass
    # calling super().save(ar) does not re-enter its own override). If the
    # exporter ever stopped binding save() on the derived type, that inherited
    # non-virtual base save is what would run, and this spelling would silently
    # checkpoint the base state only -- a partial write where the previous
    # behaviour was a loud "Unsupported type". Compare the two trees.
    def entries(archive, path="/", found=None, depth=0):
        found = [] if found is None else found
        if depth > 12:
            return found
        for child in archive.list_children(path):
            child_path = path.rstrip("/") + "/" + child
            found.append(child_path)
            if archive.is_group(child_path):
                entries(archive, child_path, found, depth + 1)
        return found

    trees = {}
    for label, write in (("explicit", lambda sim, ar: sim.save(ar)),
                         ("setitem", lambda sim, ar: ar.__setitem__("/", sim))):
        fresh = ising_c.sim(parameters)
        assert fresh.run(lambda: False)
        target = os.path.join(directory, label + ".h5")
        with hdf5.archive(target, "w") as archive:
            write(fresh, archive)
        with hdf5.archive(target, "r") as archive:
            trees[label] = sorted(entries(archive))

    assert trees["setitem"] == trees["explicit"], (
        "archive['/'] = simulation wrote a different tree than "
        "simulation.save(archive):\n"
        f"  only in explicit: {sorted(set(trees['explicit']) - set(trees['setitem']))}\n"
        f"  only in setitem:  {sorted(set(trees['setitem']) - set(trees['explicit']))}")
    # Guard against both spellings degenerating to base-only state.
    assert any(entry.endswith("/checkpoint/sweeps") for entry in trees["setitem"]), \
        trees["setitem"]

print("downstream nanobind export: ok")
