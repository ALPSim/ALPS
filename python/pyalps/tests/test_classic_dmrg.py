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


def test_warmup_configuration_places_an_unevenly_distributed_ground_state():
    # Fields alone, up spins on the first quarter of 32 sites: the uniform
    # warm-up misplaces spins the sweeps cannot move back; the configuration
    # (local state 0 up, 1 down) gives the exact wall.
    L, wall = 32, 8
    xx = models.xxz(lat.chain(L, bc="open"), Jxy=1.0, Jz=0.0)
    site = xx.site_for(0)
    fields = models.TermSum()
    for i in range(L):
        fields += (-1.0 if i < wall else 1.0) * site.Sz(i)
    op = models.compile_observable(xx, fields)
    configuration = [0] * wall + [1] * (L - wall)
    result = dmrg.solve(op, xx.physical_basis(), sector={"Sz": wall - L / 2}, sweeps=dmrg.fixed_sweeps(4, 32),
                        warmup_configuration=configuration)
    assert result.energy == pytest.approx(-0.5 * L, abs=1e-10)  # only the wall reaches -L/2
    assert list(result.options.warmup_configuration) == configuration
    with pytest.raises(ValueError, match="warmup_configuration has Sz = -8, but the sector asks for 0"):
        dmrg.solve(op, xx.physical_basis(), sector={"Sz": 0}, warmup_configuration=configuration)
    with pytest.raises(ValueError, match="one state per site"):
        dmrg.solve(op, xx.physical_basis(), sector={"Sz": 0}, warmup_configuration=[0, 1])


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
    # Operator names and located terms as measurements, on catalog and kit models alike.
    named = dmrg.solve(ring, sector={"Sz": 0.0}, measurements={"Sz": "Sz"}, **exact_sweeps)
    np.testing.assert_allclose(named.measurements["Sz"].values[0], 0.0, atol=1e-8)
    assert len(named.measurements["Sz"].values[0]) == 8
    spin = models.spin()
    kit = models.heisenberg(models.lattice.chain(8))
    located = dmrg.solve(kit, sector={"Sz": 0.0}, expectation={"SzSz": spin.Sz(0) * spin.Sz(1)}, **exact_sweeps)
    ground = ed.ground_state(kit, sector={"Sz": 0.0}, observables={"SzSz": spin.Sz(0) * spin.Sz(1)})
    assert located.energy == pytest.approx(ground.energy, abs=1e-8)
    assert located.measurements["SzSz"].values[0, 0] == pytest.approx(ground.observables["SzSz"], abs=1e-7)
    with pytest.raises(TypeError, match="needs a model"):
        dmrg.make_measurements({"Sz": "Sz"})


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


# -- keep_mps: the final states as matrix product states ----------------------------

EXACT = {"sweeps": dmrg.fixed_sweeps(4, 128), "warmup_states": 128, "lanczos_tolerance": 1e-13}


def assert_is_eigenvector(mps, reference, column):
    """The MPS on the ED sector states is the ED eigenvector up to a phase, with no weight elsewhere."""
    problem = reference.problem
    amplitudes = mps.amplitude(problem.states)
    assert np.linalg.norm(amplitudes) == pytest.approx(1.0, abs=1e-10)
    assert abs(np.vdot(reference.eigenvectors[:, column], amplitudes)) == pytest.approx(1.0, abs=1e-8)
    assert mps.norm() == pytest.approx(1.0, abs=1e-10)


def test_kept_mps_of_kit_models_are_the_exact_states():
    model = models.hubbard(lat.chain(8, bc="open"), t=1.0, U=4.0, mu=0.3)
    sector = {"N": 7.0, "Sz": 0.5}
    result = dmrg.solve(model, sector=sector, num_eigenvalues=2, keep_mps=True, **EXACT)
    reference = exact(model, count=2, **sector)
    assert len(result.mps) == 2
    first, second = result.mps
    assert isinstance(first, dmrg.MatrixProductState) and len(first) == 8
    assert first.center == 4 and first.local_dimensions == [4] * 8 and not first.is_complex
    tensors = first.tensors
    assert [t.shape[0] for t in tensors[1:]] == first.bond_dimensions
    assert tensors[0].shape[0] == 1 and tensors[-1].shape[2] == 1 and tensors[0].dtype == np.float64
    for column, mps in enumerate(result.mps):
        assert_is_eigenvector(mps, reference, column)
    assert abs(first.overlap(second)) < 1e-8
    dense = first.to_dense()
    assert dense.shape == (4 ** 8,)
    assert np.vdot(dense, dense).real == pytest.approx(1.0, abs=1e-10)
    np.testing.assert_allclose(first.amplitude(reference.problem.states[0]), first.amplitude(reference.problem.states[:1])[0])
    with pytest.raises(ValueError, match="local index"):
        first.amplitude([4] * 8)
    with pytest.raises(ValueError, match="max_dimension"):
        first.to_dense(max_dimension=100)
    assert dmrg.solve(model, sector=sector, **EXACT).mps == []


def test_kept_mps_in_complex_arithmetic_and_with_stored_products():
    site = models.spinless_fermion()
    model = models.Hamiltonian(lat.nnn_chain(8, bc="open"), site, conserved_quantum_numbers=["N"])
    hop = lambda i, j: models.located(site, "c_dag", i) * models.located(site, "c", j)  # noqa: E731
    model.add_bond(-np.exp(0.4j), hop, bond_type=0, add_hc=True)
    model.add_bond(-0.5, hop, bond_type=1, add_hc=True)
    result = dmrg.solve(model, sector={"N": 3.0}, arithmetic="complex", keep_mps=True, **EXACT)
    (mps,) = result.mps
    assert mps.is_complex and mps.tensors[0].dtype == np.complex128
    assert_is_eigenvector(mps, exact(model, N=3.0), 0)
    # Two targets with a two-site measurement store products: the final sweep closes at the right end.
    chain = models.heisenberg(lat.chain(8, bc="open"))
    stored = dmrg.solve(chain, sector={"Sz": 0.0}, num_eigenvalues=2, keep_mps=True,
                        measurements={"SzSz": dmrg.pair_terms(chain, "Sz", "Sz")}, **EXACT)
    reference = exact(chain, count=2, Sz=0.0)
    assert [m.center for m in stored.mps] == [6, 6]
    for column, kept in enumerate(stored.mps):
        assert_is_eigenvector(kept, reference, column)


def test_kept_mps_of_catalog_models_use_the_catalog_basis():
    spin = models.CatalogModel({"MODEL": "spin", "LATTICE": "open chain lattice", "L": 8, "J": 1.0,
                                "local_S": 1.0, "Sz_total": 1, "CONSERVED_QUANTUMNUMBERS": "Sz"})
    result = dmrg.solve(spin, keep_mps=True, **EXACT)
    assert result.mps[0].local_dimensions == [3] * 8
    assert_is_eigenvector(result.mps[0], ed.prepare_product(spin).eigenstates(count=1), 0)
    hubbard = models.CatalogModel({"MODEL": "fermion Hubbard", "LATTICE": "open chain lattice", "L": 6,
                                   "t": 1.0, "U": 4.0, "Nup_total": 3, "Ndown_total": 2,
                                   "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown"})
    result = dmrg.solve(hubbard, keep_mps=True, **EXACT)
    assert_is_eigenvector(result.mps[0], ed.prepare_product(hubbard).eigenstates(count=1), 0)


def test_kept_mps_of_a_truncated_run_is_its_own_state():
    model = models.heisenberg(lat.chain(16, bc="open"))
    result = dmrg.solve(model, sector={"Sz": 0.0}, sweeps=dmrg.fixed_sweeps(3, 12), warmup_states=12, keep_mps=True)
    (mps,) = result.mps
    assert mps.center == 8 and max(mps.bond_dimensions) <= 24
    assert mps.norm() == pytest.approx(1.0, abs=1e-10)
    # Its energy is the run's: <psi|H|psi> on the dense vector.
    psi = mps.to_dense()
    H = model.to_sparse_matrix()
    assert np.vdot(psi, H @ psi).real == pytest.approx(result.energy, abs=1e-8)


@pytest.mark.parametrize("n, n_alpha, n_beta, kept", [(2, 1, 1, 128), (5, 3, 2, 128), (6, 3, 3, 8)])
def test_kept_molecular_mps_is_the_ci_vector(n, n_alpha, n_beta, kept):
    h1, eri = symmetric_integrals(n, seed=5)
    result = dmrg.solve_molecular(h1, eri, n_alpha, n_beta, sweeps=dmrg.fixed_sweeps(4, kept), warmup_states=kept,
                                  lanczos_tolerance=1e-12, ci_vector=True, keep_mps=True)
    (mps,) = result.mps
    assert len(mps) == n and mps.local_dimensions == [4] * n
    alpha, beta = dmrg.ci_occupations(result)
    amplitudes = mps.amplitude(alpha + 2 * beta)  # spinful_fermion states 0, up, down, up-down
    sign = np.sign(np.vdot(amplitudes, result.ci_coefficients))
    np.testing.assert_allclose(sign * amplitudes, result.ci_coefficients, atol=1e-10)
    assert mps.norm() == pytest.approx(1.0, abs=1e-10)


# -- evolve: adaptive time-dependent DMRG ---------------------------------------------

def exact_evolution(initial, final, times, **sector):
    """The ED ground state of ``initial`` and exp(-i H_final t) applied to it, on the sector basis."""
    reference = exact(initial, **sector)
    psi0 = reference.eigenvectors[:, 0]
    energies, vectors = np.linalg.eigh(final.to_dense_matrix(sector))
    coefficients = vectors.conj().T @ psi0
    return reference, [vectors @ (np.exp(-1j * energies * t) * coefficients) for t in times]


def wall_magnetization(n, t):
    """<Sz> a distance n right of an XX domain wall on the infinite chain (Antal et al., PRE 59, 4912)."""
    from scipy.special import jv

    return -0.5 * sum(jv(i, t) ** 2 for i in range(1 - n, n))


def test_evolve_melts_the_xx_domain_wall_as_the_bessel_solution():
    L, half = 16, 8
    xx = models.xxz(lat.chain(L, bc="open"), Jxy=1.0, Jz=0.0)
    site = xx.site_for(0)
    # Up spins left, down spins right: the product ground state of the
    # fields, which the warm-up configuration places.
    wall = models.TermSum()
    for i in range(L):
        wall += (-1.0 if i < half else 1.0) * site.Sz(i)
    run = dmrg.evolve(wall, xx, time_step=0.05, steps=60, states=40, truncation_error=1e-10, sector={"Sz": 0},
                      warmup_configuration=[0] * half + [1] * half, measure_steps=5, measurements={"Sz": "Sz"})
    assert run.initial.energy == pytest.approx(-0.5 * L, abs=1e-10)
    np.testing.assert_allclose(run.steps, np.arange(0, 61, 5))
    np.testing.assert_allclose(run.times, 0.05 * run.steps)
    profiles = run.measurements["Sz"].values
    assert profiles.shape == (13, L)
    np.testing.assert_allclose(profiles[0], [0.5] * half + [-0.5] * half, atol=1e-10)
    np.testing.assert_allclose(profiles.sum(axis=1), 0.0, atol=1e-10)
    for k, t in enumerate(run.times):
        for n in (1, 2, 3):
            assert profiles[k, half + n - 1] == pytest.approx(wall_magnetization(n, t), abs=5e-4), (t, n)
    # The front has not reached the ends at t = 3; the energy is conserved.
    assert abs(profiles[-1, 0]) > 0.499 and abs(profiles[-1, -1]) > 0.499
    np.testing.assert_allclose(run.energies, run.energies[0], atol=1e-4)
    assert run.truncation_errors.shape == (61,) and run.truncation_errors.sum() < 1e-6


def test_evolve_states_and_loschmidt_echo_match_exact_evolution():
    chain = lat.chain(8, bc="open")
    initial, final = models.xxz(chain, Jxy=1.0, Jz=0.0), models.xxz(chain, Jxy=1.0, Jz=2.0, h=0.2)
    run = dmrg.evolve(initial, final, time_step=0.02, steps=60, states=64, sector={"Sz": 0.0}, measure_steps=6,
                      keep_mps=True, **EXACT)
    reference, states = exact_evolution(initial, final, run.times, Sz=0.0)
    psi0 = reference.eigenvectors[:, 0]
    assert run.initial.energy == pytest.approx(reference.eigenvalues[0], abs=1e-9)
    assert len(run.mps) == len(run.times) == 11 and run.mps[0].is_complex
    amplitudes = [mps.amplitude(reference.problem.states) for mps in run.mps]
    for psi, kept in zip(states, amplitudes):
        assert abs(np.vdot(psi, kept)) == pytest.approx(1.0, abs=1e-4)
    # <psi(0)|psi(t)> is free of the initial state's phase.
    expected = np.array([np.vdot(psi0, psi) for psi in states])
    np.testing.assert_allclose(dmrg.overlaps(run), expected, atol=2e-4)
    np.testing.assert_allclose(dmrg.loschmidt_echo(run), np.abs(expected) ** 2, atol=2e-4)
    assert dmrg.loschmidt_echo(run)[0] == pytest.approx(1.0, abs=1e-12)
    H = final.to_dense_matrix({"Sz": 0.0})
    assert run.energies[0] == pytest.approx(np.vdot(psi0, H @ psi0).real, abs=1e-9)
    np.testing.assert_allclose(run.energies, run.energies[0], atol=1e-4)
    # Overlaps with another state on the same sites.
    ground = dmrg.solve(final, sector={"Sz": 0.0}, keep_mps=True, **EXACT).mps[0]
    final_ground = exact(final, Sz=0.0).eigenvectors[:, 0]
    np.testing.assert_allclose(np.abs(dmrg.overlaps(run, ground)),
                               [abs(np.vdot(final_ground, psi)) for psi in states], atol=2e-4)


def test_evolve_trotter_error_scales_as_the_step_squared():
    chain = lat.chain(8, bc="open")
    initial, final = models.xxz(chain, Jxy=1.0, Jz=0.0), models.xxz(chain, Jxy=1.0, Jz=1.0, h=0.2)
    reference, (target,) = exact_evolution(initial, final, [1.0], Sz=0.0)
    errors = []
    for dt in (0.2, 0.1, 0.05):
        steps = round(1.0 / dt)
        run = dmrg.evolve(initial, final, time_step=dt, steps=steps, measure_steps=[steps], keep_mps=True, states=64,
                          sector={"Sz": 0.0}, **EXACT)
        psi = run.mps[0].amplitude(reference.problem.states)
        overlap = np.vdot(psi, target)
        errors.append(np.linalg.norm(target - overlap / abs(overlap) * psi))
    ratios = np.array(errors[:-1]) / np.array(errors[1:])
    assert np.all((ratios > 3.5) & (ratios < 4.5)), errors


def test_evolve_in_imaginary_time_reaches_the_ground_state():
    chain = lat.chain(8, bc="open")
    initial, final = models.xxz(chain, Jxy=1.0, Jz=0.0), models.heisenberg(chain)
    e0 = exact(final, Sz=0.0).eigenvalues[0]
    run = dmrg.evolve(initial, final, time_step=0.05, steps=200, imaginary=True, measure_steps=50, states=64,
                      sector={"Sz": 0.0}, keep_mps=True, **EXACT)
    # Down to the Trotter operator's fixed point, whose energy lies O(dt^2) above E0.
    assert np.all(np.diff(run.energies) < 1e-8)
    assert run.energies[-1] == pytest.approx(e0, abs=2e-3)
    assert run.times[-1] == pytest.approx(10.0)
    # Imaginary time keeps the options' (real) arithmetic.
    assert run.initial.options.arithmetic == dmrg.Arithmetic.real and not run.mps[-1].is_complex
    complex_run = dmrg.evolve(initial, final, time_step=0.05, steps=200, imaginary=True, measure_steps=[200],
                              states=64, sector={"Sz": 0.0}, arithmetic="complex", **EXACT)
    assert complex_run.energies[0] == pytest.approx(run.energies[-1], abs=1e-9)


def test_evolve_drives_are_sampled_at_the_step_midpoints():
    chain = lat.chain(8, bc="open")
    model = models.xxz(chain, Jxy=1.0, Jz=0.5)
    site = model.site_for(0)
    staggered = models.TermSum()
    for i in range(8):
        staggered += (-1.0) ** i * site.Sz(i)
    dt, steps = 0.05, 12

    def field(t):
        return 0.5 * np.sin(2.0 * t)

    common = dict(time_step=dt, steps=steps, states=64, sector={"Sz": 0.0}, measurements={"Sz": "Sz"},
                  keep_mps=True, **EXACT)
    by_function = dmrg.evolve(model, model, drives=[(staggered, field)], **common)
    samples = [field((n + 0.5) * dt) for n in range(steps)]
    by_values = dmrg.evolve(model, model, drives=[dmrg.Drive(models.compile_observable(model, staggered), samples)],
                            **common)
    np.testing.assert_allclose(by_function.measurements["Sz"].values, by_values.measurements["Sz"].values, atol=1e-12)
    # The ground state stays put without the drive and responds to it.
    static = dmrg.evolve(model, model, **common)
    np.testing.assert_allclose(static.measurements["Sz"].values, 0.0, atol=1e-8)
    assert np.abs(by_function.measurements["Sz"].values[-1]).max() > 1e-3
    # A piecewise-constant exact evolution with the same samples.
    reference = exact(model, Sz=0.0)
    psi = reference.eigenvectors[:, 0].astype(complex)
    H = model.to_dense_matrix({"Sz": 0.0})
    D = reference.problem.operator_matrix(staggered).toarray()
    for n in range(steps):
        energies, vectors = np.linalg.eigh(H + samples[n] * D)
        psi = vectors @ (np.exp(-1j * energies * dt) * (vectors.conj().T @ psi))
    kept = by_function.mps[-1].amplitude(reference.problem.states)
    assert abs(np.vdot(psi, kept)) == pytest.approx(1.0, abs=1e-4)
    # The energies are the final Hamiltonian's, without the drive.
    assert by_function.energies[-1] == pytest.approx(np.vdot(kept, H @ kept).real, abs=1e-9)


def test_evolve_on_catalog_models_and_compiled_operators():
    params = {"MODEL": "spin", "LATTICE": "open chain lattice", "L": 8, "J": 1.0, "Sz_total": 0,
              "CONSERVED_QUANTUMNUMBERS": "Sz"}
    initial = models.CatalogModel(params)
    halved = models.CatalogModel({**params, "J": 0.5})
    # An eigenstate of H is one of H / 2: the energy halves and the echo stays one.
    run = dmrg.evolve(initial, halved, time_step=0.1, steps=5, keep_mps=True, states=64, **EXACT)
    assert run.initial.options.sector == {"Sz": 0.0}
    np.testing.assert_allclose(run.energies, 0.5 * run.initial.energy, atol=1e-6)
    np.testing.assert_allclose(dmrg.loschmidt_echo(run), 1.0, atol=1e-6)
    # Compiled operators on an explicit basis.
    chain = lat.chain(8, bc="open")
    a, b = models.xxz(chain, Jxy=1.0, Jz=0.0), models.xxz(chain, Jxy=1.0, Jz=2.0)
    compiled = dmrg.evolve(models.compile_hamiltonian(a).op, models.compile_hamiltonian(b).op, a.physical_basis(),
                           time_step=0.05, steps=4, states=64, sector={"Sz": 0.0}, **EXACT)
    direct = dmrg.evolve(a, b, time_step=0.05, steps=4, states=64, sector={"Sz": 0.0}, **EXACT)
    np.testing.assert_allclose(compiled.energies, direct.energies, atol=1e-10)
    evolution = dmrg.make_evolution(time_step=0.05, steps=4, states=64)
    assert dmrg.evolve(a, b, evolution=evolution, sector={"Sz": 0.0}, **EXACT).energies == pytest.approx(direct.energies)


def test_evolve_rejects_invalid_requests():
    chain = lat.chain(8, bc="open")
    model = models.xxz(chain, Jxy=1.0, Jz=0.0)
    with pytest.raises(TypeError, match="time_step"):
        dmrg.evolve(model, model, steps=3)
    with pytest.raises(TypeError, match="unknown evolution field"):
        dmrg.make_evolution(time_step=0.1, stpes=3)
    with pytest.raises(TypeError, match="basis"):
        dmrg.evolve(models.compile_hamiltonian(model).op, models.compile_hamiltonian(model).op, time_step=0.1, steps=1)
    with pytest.raises(ValueError, match="nearest-neighbour"):
        dmrg.evolve(model, models.heisenberg(lat.nnn_chain(8, bc="open")), time_step=0.1, steps=1, sector={"Sz": 0.0})
    with pytest.raises(ValueError, match="measure_steps"):
        dmrg.evolve(model, model, time_step=0.1, steps=2, measure_steps=0)
    with pytest.raises(ValueError, match="after the last step"):
        dmrg.evolve(model, model, time_step=0.1, steps=2, measure_steps=[3], sector={"Sz": 0.0})
    with pytest.raises(ValueError):
        dmrg.evolve(model, models.heisenberg(lat.chain(8, bc="open"), S=1.0), time_step=0.1, steps=1)
    run = dmrg.evolve(model, model, time_step=0.1, steps=2, measure_steps=[1, 2], sector={"Sz": 0.0}, keep_mps=True)
    with pytest.raises(ValueError, match="measure step 0"):
        dmrg.loschmidt_echo(run)
    assert dmrg.loschmidt_echo(run, run.mps[0]) == pytest.approx([1.0, 1.0], abs=1e-6)
    with pytest.raises(ValueError, match="keep_mps"):
        dmrg.overlaps(dmrg.evolve(model, model, time_step=0.1, steps=1, sector={"Sz": 0.0}))

    class Stop(Exception):
        pass

    def abort(iteration):
        if np.isnan(iteration.energy):  # the evolution's iterations carry no energy
            raise Stop

    with pytest.raises(Stop):
        dmrg.evolve(model, model, time_step=0.1, steps=3, sector={"Sz": 0.0}, callback=abort)
