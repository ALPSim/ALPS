"""Mapping mutation must not expose invalid C++ map nodes to Python."""

import os
import subprocess
import sys

import pytest


_SETUP = """
from pyalps import ngs

class Simulation(ngs.mcbase):
    def update(self): pass
    def measure(self): pass
    def fraction_completed(self): return 1.0

sim = Simulation({'SEED': 1})
sim.measurements.createRealObservable('x')
sim.measurements['x'] << 1.0 << 3.0
mapping = sim.measurements if kind == 'observables' else ngs.collectResults(sim)

def mean(value):
    return (ngs.observable2result(value) if kind == 'observables' else value).mean
"""


def _run(kind, code):
    # A regression is a native use-after-free. Isolate it so a crash reports
    # an ordinary test failure instead of taking down the whole suite.
    completed = subprocess.run(
        [sys.executable, "-c", f"kind = {kind!r}\n" + _SETUP + code],
        env={**os.environ, "MallocScribble": "1"},
        capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.parametrize("kind", ["observables", "results"])
@pytest.mark.parametrize("remove", ["del mapping['x']", "mapping.clear()", "mapping.pop('x')"])
def test_values_survive_removal(kind, remove):
    _run(kind, f"""
value = mapping['x']
assert mean(value) == 2.0
{remove}
assert mean(value) == 2.0
del mapping, sim
import gc
gc.collect()
assert mean(value) == 2.0
""")


@pytest.mark.parametrize("kind", ["observables", "results"])
def test_pop_returns_usable_value(kind):
    _run(kind, """
value = mapping.pop('x')
assert mean(value) == 2.0
assert not mapping
""")


@pytest.mark.parametrize("kind", ["observables", "results"])
def test_assignment_and_update_replace_existing_values(kind):
    _run(kind, """
replacement = ngs.createRealObservable('replacement')
replacement << 7.0 << 9.0
if kind == 'results':
    replacement = ngs.observable2result(replacement)
mapping['x'] = replacement
assert mean(mapping['x']) == 8.0
mapping.update({'x': replacement})
assert mean(mapping['x']) == 8.0
mapping['new'] = replacement
assert mean(mapping['new']) == 8.0
""")


def test_observable_samples_still_update_the_container():
    _run('observables', """
value = mapping['x']
value << 8.0
assert mean(mapping['x']) == 4.0
""")


def test_observable_merge_still_updates_the_container():
    _run('observables', """
other = ngs.createRealObservable('x')
other << 7.0 << 9.0
mapping['x'].merge(other)
assert mean(mapping['x']) == 5.0
""")


@pytest.mark.parametrize("kind", ["observables", "results"])
def test_mapping_iterator_survives_mutation(kind):
    _run(kind, """
iterator = iter(mapping)
mapping.clear()
del mapping, sim
assert list(iterator) == ['x']
""")


def test_params_iterator_survives_mutation():
    _run('observables', """
parameters = ngs.params({'first': 1, 'second': 2})
iterator = iter(parameters)
parameters.clear()
parameters['new'] = 3
del parameters
assert list(iterator) == ['first', 'second']
""")
