# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""pyalps.classic_dmrg against pyalps.ed: energies, measurements, targets,
iteration observation and the option surface."""

import itertools

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


def catalog_reference(op, basis, sector, count):
    problem = ed.prepare_projected(op, ed.projected_spaces(basis), sector=sector)
    return problem.eigenstates(count=count).eigenvalues


def test_catalog_models_solve_like_exact_diagonalization():
    # The operators the XML dmrg application solves. The catalog's fermion
    # Hubbard basis has Nup and Ndown but no N, so both carry the fermion
    # parity; the periodic chain's wrap bond arrives in descending site order.
    exact_sweeps = {"sweeps": dmrg.fixed_sweeps(4, 200), "warmup_states": 200, "lanczos_tolerance": 1e-13}
    hubbard = models.CatalogModel({"MODEL": "fermion Hubbard", "LATTICE": "open chain lattice", "L": 6,
                                   "t": 1.0, "U": 4.0})
    op = hubbard.compile_hamiltonian()
    basis = hubbard.physical_basis(conserved_quantum_numbers=["Nup", "Ndown"])
    assert basis.conserved_quantum_numbers == ["Ndown", "Nup"]
    assert hubbard.physical_basis().conserved_quantum_numbers == hubbard.conserved_quantum_numbers()
    sector = {"Nup": 3.0, "Ndown": 3.0}
    result = dmrg.solve(op, basis, sector=sector, num_eigenvalues=2, **exact_sweeps)
    np.testing.assert_allclose(result.energies, catalog_reference(op, basis, sector, 2), atol=1e-8)

    ring = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 8, "J": 1.0, "local_S": 0.5,
                                "CONSERVED_QUANTUMNUMBERS": "Sz"})
    op = ring.compile_hamiltonian()
    basis = ring.physical_basis()
    result = dmrg.solve(op, basis, sector={"Sz": 0.0}, **exact_sweeps)
    assert result.energy == pytest.approx(catalog_reference(op, basis, {"Sz": 0.0}, 1)[0], abs=1e-8)

    # The catalog itself: its compiled Hamiltonian and basis, in its Sz_total sector unless told otherwise.
    ring.set_parameters({"Sz_total": 1})
    direct = dmrg.solve(ring, **exact_sweeps)
    assert direct.options.sector == {"Sz": 1.0}
    assert direct.energy == pytest.approx(catalog_reference(op, basis, {"Sz": 1.0}, 1)[0], abs=1e-8)
    assert dmrg.solve(ring, sector={"Sz": 0.0}, **exact_sweeps).energy == pytest.approx(result.energy, abs=1e-8)
    with pytest.raises(TypeError):
        dmrg.solve(ring, basis)


def symmetric_integrals(n, seed):
    rng = np.random.default_rng(seed)
    h1 = rng.normal(size=(n, n))
    h1 = 0.5 * (h1 + h1.T)
    eri = 0.3 * rng.normal(size=(n, n, n, n))
    eri = (eri + eri.transpose(1, 0, 2, 3)) / 2
    eri = (eri + eri.transpose(0, 1, 3, 2)) / 2
    eri = (eri + eri.transpose(2, 3, 0, 1)) / 2
    return h1, eri


def fock_reference(h1, eri, n_alpha, n_beta):
    """Full CI on the Jordan-Wigner Fock space: energy and spin-resolved RDMs."""
    energy, _, expect, n = fock_ground_state(h1, eri, n_alpha, n_beta)
    one = [np.array([[expect(("+", q, s), ("-", p, s)) for q in range(n)] for p in range(n)]) for s in (0, 1)]
    two = np.zeros((n, n, n, n))
    for p, q, r, t in itertools.product(range(n), repeat=4):
        two[p, q, r, t] = sum(expect(("+", p, s1), ("+", r, s2), ("-", t, s2), ("-", q, s1))
                              for s1, s2 in itertools.product((0, 1), repeat=2))
    return energy, one[0], one[1], two


def fock_ground_state(h1, eri, n_alpha, n_beta):
    """The ground state on the Jordan-Wigner Fock space of modes 2p + spin (spin 0 = alpha).

    Returns the energy, psi over all 2**(2n) occupation states (mode 0 the
    most significant bit), expect(*ladder) for ladder factors ("+"/"-", p,
    spin) applied right to left, and n.
    """
    from scipy import sparse

    n = h1.shape[0]
    modes = 2 * n  # mode 2p + spin, spin 0 = alpha; mode 0 is the leading tensor factor
    z, a = sparse.diags([1.0, -1.0]), sparse.csr_matrix([[0.0, 1.0], [0.0, 0.0]])
    def annihilate(m):
        out = sparse.identity(1, format="csr")
        for k in range(modes):
            out = sparse.kron(out, z if k < m else a if k == m else sparse.identity(2), format="csr")
        return out
    c = [annihilate(m) for m in range(modes)]
    cd = [x.T.tocsr() for x in c]
    occupied = lambda i, m: (i >> (modes - 1 - m)) & 1
    keep = [i for i in range(2 ** modes)
            if sum(occupied(i, 2 * p) for p in range(n)) == n_alpha
            and sum(occupied(i, 2 * p + 1) for p in range(n)) == n_beta]
    H = sparse.csr_matrix((2 ** modes, 2 ** modes))
    for p, q in itertools.product(range(n), repeat=2):
        for s in (0, 1):
            H = H + h1[p, q] * (cd[2 * p + s] @ c[2 * q + s])
    for p, q, r, t in itertools.product(range(n), repeat=4):
        for s1, s2 in itertools.product((0, 1), repeat=2):
            H = H + 0.5 * eri[p, q, r, t] * (cd[2 * p + s1] @ cd[2 * r + s2] @ c[2 * t + s2] @ c[2 * q + s1])
    values, vectors = np.linalg.eigh(H.toarray()[np.ix_(keep, keep)])
    psi = np.zeros(2 ** modes)
    psi[keep] = vectors[:, 0]
    def expect(*ops):  # <psi| ops[0] ops[1] ... |psi>, applied right to left
        phi = psi
        for kind, p, s in reversed(ops):
            phi = (cd if kind == "+" else c)[2 * p + s] @ phi
        return psi @ phi
    return values[0], psi, expect, n


@pytest.mark.parametrize("n, n_alpha, n_beta", [(2, 1, 1), (3, 2, 1), (4, 2, 2)])
def test_solve_molecular_energy_and_rdms_match_full_ci(n, n_alpha, n_beta):
    h1, eri = symmetric_integrals(n, seed=n)
    energy, alpha, beta, two = fock_reference(h1, eri, n_alpha, n_beta)
    seen = []
    result = dmrg.solve_molecular(h1, eri, n_alpha, n_beta, core_energy=0.25, sweeps=dmrg.fixed_sweeps(4, 128),
                                  warmup_states=128, lanczos_tolerance=1e-12, callback=seen.append)
    assert isinstance(result, dmrg.MolecularResult) and result.orbitals == n
    assert result.energy == pytest.approx(energy + 0.25, abs=1e-9)
    np.testing.assert_allclose(result.one_rdm_alpha, alpha, atol=1e-8)
    np.testing.assert_allclose(result.one_rdm_beta, beta, atol=1e-8)
    np.testing.assert_allclose(result.one_rdm, alpha + beta, atol=1e-8)
    assert result.two_rdm.shape == (n,) * 4
    np.testing.assert_allclose(result.two_rdm, two, atol=1e-8)
    # The RDMs reproduce the energy (PySCF make_rdm12 contraction).
    rdm_energy = 0.25 + np.einsum("pq,qp", h1, result.one_rdm) + 0.5 * np.einsum("pqrs,pqrs", eri, result.two_rdm)
    assert rdm_energy == pytest.approx(result.energy, abs=1e-9)
    assert len(seen) == len(result.iterations) and result.iteration_energies.shape == (len(seen),)


def test_solve_molecular_spin_blocks_observables_and_ci_vector_match_full_ci():
    n, n_alpha, n_beta = 4, 2, 1
    h1, eri = symmetric_integrals(n, seed=11)
    _, psi, expect, _ = fock_ground_state(h1, eri, n_alpha, n_beta)
    fs = models.spinful_fermion()
    def at(op, i):
        return models.located(fs, op, i)
    observables = {"doublon pair": at("n_updn", 0) * at("n_updn", 2), "spin flip": at("Sp", 0) * at("Sm", 3),
                   "four sites": at("n_up", 0) * at("c_dag_dn", 1) * at("n", 2) * at("c_dn", 3)}
    result = dmrg.solve_molecular(h1, eri, n_alpha, n_beta, sweeps=dmrg.fixed_sweeps(4, 128), warmup_states=128,
                                  lanczos_tolerance=1e-12, observables=observables, ci_vector=True)
    blocks = {"aa": (0, 0), "ab": (0, 1), "bb": (1, 1)}
    for name, (s1, s2) in blocks.items():
        expected = np.zeros((n,) * 4)
        for p, q, r, t in itertools.product(range(n), repeat=4):
            expected[p, q, r, t] = expect(("+", p, s1), ("+", r, s2), ("-", t, s2), ("-", q, s1))
        np.testing.assert_allclose(getattr(result, "two_rdm_" + name), expected, atol=1e-8)
    np.testing.assert_allclose(result.two_rdm_aa + result.two_rdm_ab + result.two_rdm_ab.transpose(2, 3, 0, 1)
                               + result.two_rdm_bb, result.two_rdm, atol=1e-12)
    number = lambda p, s: (("+", p, s), ("-", p, s))
    assert result.observables["doublon pair"] == pytest.approx(
        expect(*number(0, 0), *number(0, 1), *number(2, 0), *number(2, 1)), abs=1e-8)
    assert result.observables["spin flip"] == pytest.approx(
        expect(("+", 0, 0), ("-", 0, 1), ("+", 3, 1), ("-", 3, 0)), abs=1e-8)
    assert result.observables["four sites"] == pytest.approx(
        sum(expect(*number(0, 0), ("+", 1, 1), *number(2, s), ("-", 3, 1)) for s in (0, 1)), abs=1e-8)
    # CI: the reference orders an orbital's alpha mode before its beta mode,
    # ALPS the down (beta) mode first, so a doubly occupied orbital flips the
    # sign; the global sign makes the largest amplitude positive.
    alpha, beta = dmrg.ci_occupations(result)
    assert alpha.shape == beta.shape == (len(result.ci_coefficients), n) == (6 * 4, n)
    assert (alpha.sum(axis=1) == n_alpha).all() and (beta.sum(axis=1) == n_beta).all()
    modes = 2 * n
    index = [sum(int(alpha[k, p]) << (modes - 1 - 2 * p) | int(beta[k, p]) << (modes - 2 - 2 * p) for p in range(n))
             for k in range(len(alpha))]
    reference = psi[index] * (-1.0) ** (alpha & beta).sum(axis=1)
    reference *= np.sign(reference[np.argmax(np.abs(reference))])
    np.testing.assert_allclose(result.ci_coefficients, reference, atol=1e-8)
    assert np.linalg.norm(result.ci_coefficients) == pytest.approx(1.0, abs=1e-12)


def test_molecular_correlations_from_rdms_match_measured_observables():
    n = 4
    h1, eri = symmetric_integrals(n, seed=5)
    fs = models.spinful_fermion()
    names = {"n": "n", "n_up": "n_up", "n_dn": "n_dn", "Sz": "Sz"}
    observables = {}
    for i, j in itertools.product(range(n), repeat=2):
        for a, b in (("n", "Sz"), ("Sz", "Sz"), ("n_up", "n_dn")):
            observables[f"{a} {b} {i} {j}"] = models.located(fs, names[a], i) * models.located(fs, names[b], j)
        observables[f"Sp {i} {j}"] = models.located(fs, "Sm", i) * models.located(fs, "Sp", j)
        observables[f"Sm {i} {j}"] = models.located(fs, "Sp", i) * models.located(fs, "Sm", j)
        observables[f"c_up {i} {j}"] = models.located(fs, "c_dag_up", i) * models.located(fs, "c_up", j)
    result = dmrg.solve_molecular(h1, eri, 2, 2, sweeps=dmrg.fixed_sweeps(4, 64), warmup_states=64,
                                  lanczos_tolerance=1e-12, observables=observables)
    for a, b in (("n", "Sz"), ("Sz", "Sz"), ("n_up", "n_dn")):
        measured = np.array([[result.observables[f"{a} {b} {i} {j}"] for j in range(n)] for i in range(n)])
        np.testing.assert_allclose(dmrg.correlation_matrix(result, a, b), measured, atol=1e-10)
    for op in ("Sp", "Sm", "c_up"):
        measured = np.array([[result.observables[f"{op} {i} {j}"] for j in range(n)] for i in range(n)])
        np.testing.assert_allclose(dmrg.correlation_matrix(result, op), measured, atol=1e-10)
    np.testing.assert_allclose(dmrg.local_expectations(result, "n"), np.diag(result.one_rdm), atol=1e-12)
    np.testing.assert_allclose(dmrg.local_expectations(result, "n_updn"),
                               np.diag(dmrg.correlation_matrix(result, "n_up", "n_dn")), atol=1e-12)
    # A singlet: <S^2> = 0.
    assert dmrg.spin_spin_correlations(result).sum() == pytest.approx(0.0, abs=1e-8)
    with pytest.raises(ValueError, match="observables"):
        dmrg.correlation_matrix(result, "c_up", "n")
    with pytest.raises(ValueError, match="observables"):
        dmrg.local_expectations(result, "Sx")


def test_solve_molecular_rejects_invalid_observables():
    h1, eri = symmetric_integrals(5, seed=3)
    fs = models.spinful_fermion()
    for observable in (models.located(fs, "n", 5),
                       models.located(fs, "n", 0) * models.located(fs, "n", 1) * models.located(fs, "n", 2)
                       * models.located(fs, "n", 3) * models.located(fs, "n", 4),
                       models.located(fs, "Sy", 0),
                       models.located(models.spin(0.5), "Sz", 0)):
        with pytest.raises(ValueError):
            dmrg.solve_molecular(h1, eri, 2, 2, sweeps=dmrg.fixed_sweeps(1, 16), observables={"bad": observable})


def test_solve_molecular_rejects_invalid_requests():
    h1, eri = symmetric_integrals(3, seed=7)
    with pytest.raises(ValueError):
        dmrg.solve_molecular(h1, eri[0], 1, 1)  # rank
    with pytest.raises(TypeError):
        dmrg.solve_molecular(h1.astype(complex), eri, 1, 1)
    broken = eri.copy()
    broken[0, 1, 2, 2] += 1e-3
    with pytest.raises(ValueError):
        dmrg.solve_molecular(h1, broken, 1, 1)  # symmetry
    with pytest.raises(ValueError):
        dmrg.solve_molecular(h1, eri, 1, 1, sector={"N": 2.0})
    with pytest.raises(ValueError):
        dmrg.solve_molecular(h1, eri, 4, 1)  # more electrons than orbitals
    with pytest.raises(TypeError):
        dmrg.solve_molecular(h1, eri, 1, 1, callback=42)
