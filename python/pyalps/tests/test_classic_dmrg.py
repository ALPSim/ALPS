# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""pyalps.classic_dmrg against pyalps.ed: energies, measurements, targets,
iteration observation and the option surface."""

import numpy as np
import pytest

models = pytest.importorskip("pyalps.models")
ed = pytest.importorskip("pyalps.ed")
dmrg = pytest.importorskip("pyalps.classic_dmrg")
lat = models.lattice


def exact(model, count=1, **sector):
    return ed.prepare_product(model, sector=sector).eigenstates(count=count)


def test_open_chain_energy_and_measurements_match_exact_diagonalization():
    model = models.heisenberg(lat.chain(8, bc="open"))
    spin = models.spin()
    reference = exact(model, Sz=0.0)
    result = dmrg.solve(model, sweeps=dmrg.fixed_sweeps(3, 60), sector={"Sz": 0.0},
                        measurements={"Sz": dmrg.site_terms(model, "Sz"), "SzSz": dmrg.pair_terms(model, "Sz", "Sz")},
                        expectation={"H": models.compile_hamiltonian(model).op})
    assert result.energies.shape == (1,) and result.energy == pytest.approx(reference.eigenvalues[0], abs=1e-9)
    assert result.final_truncation_error >= -1e-12 and len(result.iterations) > 0  # 1 - sum(weights): rounding
    assert result.options.sector == {"Sz": 0.0} and result.options.num_eigenvalues == 1
    assert [s.right_to_left for s in result.options.sweeps] == [60, 60, 60]
    assert set(it.direction for it in result.iterations) == {"left_to_right", "right_to_left"}
    assert result.iteration_energies.shape == (len(result.iterations),)
    assert result.iteration_energies[-1] == pytest.approx(result.energy, abs=1e-9)

    local = result.measurements["Sz"]
    assert local.kind == dmrg.MeasurementKind.terms and local.sites == [[i] for i in range(8)]
    np.testing.assert_allclose(local.values, 0.0, atol=1e-8)
    pairs = result.measurements["SzSz"]
    assert pairs.values.shape == (1, 64) and pairs.sites[3 * 8 + 4] == [3, 4]
    correlation = model.operator_decomposition(models.located(spin, "Sz", 3) * models.located(spin, "Sz", 4)).op
    assert pairs.values[0, 3 * 8 + 4] == pytest.approx(reference.expectation_values(correlation)[0].real, abs=1e-8)
    on_site = model.operator_decomposition(models.located(spin, "Sz", 2) * models.located(spin, "Sz", 2)).op
    assert pairs.values[0, 2 * 8 + 2] == pytest.approx(reference.expectation_values(on_site)[0].real, abs=1e-8)
    energy = result.measurements["H"]
    assert energy.kind == dmrg.MeasurementKind.expectation and energy.sites == [[]]
    assert energy.values[0, 0] == pytest.approx(result.energy, abs=1e-8)


def test_two_targets_in_a_fermionic_sector_with_a_callback():
    model = models.hubbard(lat.chain(6, bc="open"), t=1.0, U=4.0)
    reference = exact(model, count=2, N=6.0, Sz=0.0)
    seen = []
    result = dmrg.solve(model, sweeps=dmrg.ramped_sweeps(4, 160), sector={"N": 6.0, "Sz": 0.0}, num_eigenvalues=2,
                        measurements={"n": dmrg.site_terms(model, "n")}, callback=seen.append)
    np.testing.assert_allclose(result.energies, reference.eigenvalues, atol=1e-7)
    assert len(seen) == len(result.iterations) and all(isinstance(it, dmrg.Iteration) for it in seen)
    assert [it.index for it in seen] == [it.index for it in result.iterations]
    density = result.measurements["n"].values
    assert density.shape == (2, 6)
    np.testing.assert_allclose(density.sum(axis=1), 6.0, atol=1e-7)
    densities = model.operator_decomposition(models.located(models.spinful_fermion(), "n", 0)).op
    np.testing.assert_allclose(density[:, 0], reference.expectation_values(densities).real, atol=1e-6)
    assert result.options.sweeps == [dmrg.SweepStates(20, 40), dmrg.SweepStates(60, 80), dmrg.SweepStates(100, 120),
                                     dmrg.SweepStates(140, 160)]


def test_compiled_operator_entry_point_and_explicit_options():
    model = models.heisenberg(lat.chain(6, bc="open"))
    reference = exact(model, Sz=0.0).eigenvalues[0]
    options = dmrg.Options()
    options.sweeps = [(40, 40), (40, 40)]
    options.sector = {"Sz": 0.0}
    options.warmup_states = 16
    assert options.sweeps == [dmrg.SweepStates(40, 40)] * 2 and options.arithmetic == dmrg.Arithmetic.real
    options.arithmetic = "complex"
    assert options.arithmetic == dmrg.Arithmetic.complex
    copied = dmrg.make_options(options, arithmetic="real", lanczos_tolerance=1e-9)
    assert copied.arithmetic == dmrg.Arithmetic.real and copied.lanczos_tolerance == 1e-9 and copied.warmup_states == 16
    with pytest.raises(TypeError):
        dmrg.make_options(bogus=1)
    compiled = models.compile_hamiltonian(model).op
    basis = models.physical_basis(model)
    result = dmrg.solve(compiled, basis, options=copied, measurements=dmrg.Measurements())
    assert result.energy == pytest.approx(reference, abs=1e-9)
    assert result.options.lanczos_tolerance == pytest.approx(1e-9)
    assert dmrg.solve(model, options=options).energy == pytest.approx(reference, abs=1e-9)  # complex arithmetic


def test_abort_and_invalid_requests():
    model = models.heisenberg(lat.chain(6, bc="open"))

    class Stop(Exception):
        pass

    def abort(iteration):
        if iteration.index >= 2:
            raise Stop(iteration.index)

    with pytest.raises(Stop):
        dmrg.solve(model, sweeps=dmrg.fixed_sweeps(2, 20), callback=abort)
    # The engine state is restored: the next solve runs normally.
    assert dmrg.solve(model, sweeps=dmrg.fixed_sweeps(2, 40)).energy == pytest.approx(exact(model).eigenvalues[0], abs=1e-9)
    with pytest.raises(ValueError):
        dmrg.solve(model, sector={"Sz": 7.0})
    with pytest.raises(ValueError):
        dmrg.solve(model, sector={"N": 0.0})
    with pytest.raises(TypeError):
        dmrg.solve(model, callback=42)
    with pytest.raises(ValueError):
        dmrg.solve(models.heisenberg(lat.chain(3, bc="open")))
    assert dmrg.ramped_sweeps(2, 8) == [dmrg.SweepStates(2, 4), dmrg.SweepStates(6, 8)]
    with pytest.raises(ValueError):
        dmrg.ramped_sweeps(4, 1)
