# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""pyalps.ed: product and projected problems, the solvers, symmetry spectra,
sparse assembly and thermodynamics, cross-checked against each other and
against NumPy."""

import numpy as np
import pytest

models = pytest.importorskip("pyalps.models")
ed = pytest.importorskip("pyalps.ed")
lat = models.lattice

CHAIN_8_GROUND_ENERGY = -3.6510934089371783  # periodic S=1/2 Heisenberg ring, J = 1


def test_product_problem_on_a_charge_sector():
    model = models.heisenberg(lat.chain(8))
    problem = ed.prepare_product(model, sector={"Sz": 0.0})
    assert problem.dimension == 70 and problem.site_dims == [2] * 8
    assert problem.sector == {"Sz": 0.0} and problem.physical_basis.conserved_quantum_numbers == ["Sz"]
    states = problem.states
    assert states.shape == (70, 8) and states.dtype == np.int32
    np.testing.assert_array_equal(states.sum(axis=1), 4)
    matrix = problem.hamiltonian
    assert isinstance(matrix, ed.RealSparseMatrix) and matrix.shape == (70, 70)
    dense = matrix.toarray()
    np.testing.assert_allclose(dense, dense.T, atol=1e-14)
    assert matrix.to_scipy().nnz == matrix.nnz

    result = problem.eigenstates(count=3)
    assert result.method == "dense" and not result.sparse_requested and len(result) == 3
    assert result.eigenvalues[0] == pytest.approx(CHAIN_8_GROUND_ENERGY, abs=1e-10)
    vectors = result.eigenvectors
    assert vectors.shape == (70, 3) and vectors.dtype == np.float64
    np.testing.assert_allclose(vectors.T @ vectors, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(dense @ vectors[:, 0], result.eigenvalues[0] * vectors[:, 0], atol=1e-10)
    hamiltonian = models.compile_hamiltonian(model).op
    np.testing.assert_allclose(result.expectation_values(hamiltonian), result.eigenvalues, atol=1e-10)
    np.testing.assert_allclose(ed.expectation_values(result, hamiltonian), result.eigenvalues, atol=1e-10)
    spin = models.spin()
    correlation = model.operator_decomposition(models.located(spin, "Sz", 0) * models.located(spin, "Sz", 1)).op
    # SU(2) symmetry: <Sz Sz> on a bond is one third of the bond energy E0/N.
    assert result.expectation_values(correlation)[0].real == pytest.approx(CHAIN_8_GROUND_ENERGY / 8 / 3, abs=1e-10)
    assert problem.operator_matrix(correlation).shape == (70, 70)
    full = ed.prepare_product(model)
    assert full.dimension == 256 and full.eigenstates(count=1).eigenvalues[0] == pytest.approx(CHAIN_8_GROUND_ENERGY, abs=1e-10)
    with pytest.raises(ValueError):
        ed.prepare_product(model, sector={"N": 0.0})


def test_lanczos_dense_and_scipy_paths_agree():
    model = models.heisenberg(lat.chain(12))
    problem = ed.prepare_product(model, sector={"Sz": 0.0})
    assert problem.dimension == 924
    dense = problem.eigenstates(count=3, method="dense")
    lanczos = problem.eigenstates(count=3, method="lanczos")
    assert dense.method == "dense" and lanczos.method == "lanczos" and lanczos.sparse_requested
    np.testing.assert_allclose(lanczos.eigenvalues, dense.eigenvalues, atol=1e-8)
    assert lanczos.iterations > 0 and lanczos.tolerance == pytest.approx(1e-10)
    assert lanczos.residual_norms.shape == (3,) and lanczos.residual_norms.max() < 1e-7
    automatic = problem.eigenstates(count=3)
    assert automatic.method == "dense"  # below the automatic Lanczos threshold

    matrix = problem.hamiltonian.to_scipy()
    direct = ed.lanczos(matrix, count=3)
    np.testing.assert_allclose(direct.eigenvalues, dense.eigenvalues, atol=1e-8)
    assert direct.eigenvectors.shape == (924, 3) and direct.nconv == 3 and direct.runs >= 1
    values_only = ed.lanczos_csr(matrix.indptr, matrix.indices, matrix.data, 924, 2, want_vectors=False)
    assert values_only.eigenvectors is None
    np.testing.assert_allclose(values_only.eigenvalues, dense.eigenvalues[:2], atol=1e-8)
    with pytest.raises(ValueError):
        ed.lanczos(matrix, count=1, which="SM")

    array = matrix.toarray()
    solved = ed.dense_eigensystem(array, count=3)
    np.testing.assert_allclose(solved.eigenvalues, dense.eigenvalues, atol=1e-10)
    np.testing.assert_allclose(array @ solved.eigenvectors, solved.eigenvectors * solved.eigenvalues, atol=1e-9)
    assert ed.dense_eigensystem(array, eigenvectors=False).eigenvectors is None
    assert len(ed.dense_eigensystem(array, eigenvectors=False)) == 924
    complex_array = array.astype(complex)
    np.testing.assert_allclose(ed.dense_eigensystem(complex_array, count=2).eigenvalues, dense.eigenvalues[:2], atol=1e-10)
    assert ed.select_eigenvalues([3.0, -1.0, 2.0], 2, "LA") == [2, 0]
    with pytest.raises(ValueError):
        ed.validate_eigensolver_controls("XX")


def test_momentum_projection_and_model_spectrum():
    model = models.heisenberg(lat.chain(8))
    assert ed.select_momentum([8], [np.pi]).index == 4
    ground = ed.prepare_projected(model, sector={"Sz": 0.0}, momentum=[0.0])
    assert ground.charges == ["Sz"] and 0 < ground.dimension < 70
    assert ground.basis.momentum == [0.0] and ground.basis.parent_dimension == 70
    assert ground.hamiltonian.shape == (ground.dimension, ground.dimension)
    result = ground.eigenstates(count=2)
    assert result.eigenvalues[0] == pytest.approx(CHAIN_8_GROUND_ENERGY, abs=1e-10)
    assert result.eigenvectors.shape == (ground.dimension, 2) and result.eigenvectors.dtype == np.complex128
    assert result.spectrum.momentum == [0.0] and result.spectrum.block_dimension == ground.dimension
    hamiltonian = models.compile_hamiltonian(model).op
    np.testing.assert_allclose(result.expectation_values(hamiltonian), result.eigenvalues, atol=1e-10)
    expansion = ed.expand_projected_state(ground.basis, result.eigenvectors[:, 0])
    assert expansion.basis.dimension == 70 and np.linalg.norm(expansion.amplitudes) == pytest.approx(1.0)
    at_pi = ed.prepare_projected(model, sector={"Sz": 0.0}, momentum=[np.pi])
    assert at_pi.eigenstates(count=1).eigenvalues[0] > result.eigenvalues[0]

    spectrum = ed.model_spectrum(model, all_momenta=True)
    assert spectrum.charges == ["Sz"] and spectrum.directions == [0]
    assert sum(len(block) for block in spectrum.blocks) == 256
    assert sum(block.block_dimension for block in spectrum.blocks) == 256
    assert {tuple(block.momentum) for block in spectrum.blocks} <= {(2 * np.pi * j / 8,) for j in range(8)}
    reference = ed.prepare_product(model).eigenstates(count=256, method="dense").eigenvalues
    np.testing.assert_allclose(spectrum.eigenvalues, reference, atol=1e-9)
    truncated = ed.model_spectrum(model, sector={"Sz": 0.0}, momentum=[np.pi], count=2)
    assert len(truncated.blocks) == 1 and len(truncated.blocks[0]) == 2
    assert truncated.sector == {"Sz": 0.0} and truncated.momentum == [np.pi]
    charge_only = ed.model_spectrum(model, sector={"Sz": 1.0})
    assert len(charge_only.blocks) == 1 and charge_only.blocks[0].block_dimension == 56


def test_explicit_projected_bases_and_operators():
    model = models.heisenberg(lat.chain(6))
    basis = models.physical_basis(model)
    spaces = ed.projected_spaces(basis)
    assert spaces.charges == ["Sz"] and spaces.num_sites == 6
    assert spaces.local_spaces[0].parity == [1.0, 1.0] and spaces.local_spaces[0].quantum_numbers == [[0.5, -0.5]]
    translation = ed.SpectrumTranslation(model.lattice.translation_generators()[0], 6)
    # Six sites: N/2 is odd, so the ground state carries momentum pi (index 3).
    source = ed.make_projected_basis(spaces.local_spaces, [translation], charge_targets=[0.0], momentum_index=3)
    target = ed.make_projected_basis(spaces.local_spaces, [translation], charge_targets=[1.0], momentum_index=3)
    assert source.num_sites == 6 and source.parent_dimension == 20 and source.representatives.shape == (source.dimension, 6)
    assert source.charge_targets == [0.0] and source.finite_group is None and len(source.orbit_sizes) == source.dimension
    hamiltonian = models.compile_hamiltonian(model).op
    problem = ed.prepare_projected(hamiltonian, source)
    assert problem.charges == [] and problem.dimension == source.dimension
    energies = problem.eigenstates(count=2).eigenvalues
    assert energies[0] == pytest.approx(ed.prepare_product(model).eigenstates(count=1).eigenvalues[0], abs=1e-10)
    via_spaces = ed.prepare_projected(hamiltonian, spaces, [translation], sector={"Sz": 0.0}, momentum=[np.pi])
    np.testing.assert_allclose(via_spaces.eigenstates(count=2).eigenvalues, energies, atol=1e-12)
    reduced = ed.projected_sparse_operator(hamiltonian, source, source)
    np.testing.assert_allclose(reduced.toarray(), problem.hamiltonian.toarray(), atol=1e-12)
    np.testing.assert_allclose(reduced.toarray(), reduced.toarray().conj().T, atol=1e-12)
    state = problem.eigenstates(count=1).eigenvectors[:, 0]
    np.testing.assert_allclose(ed.apply_projected_operator(hamiltonian, source, source, state), energies[0] * state, atol=1e-10)
    assert ed.expect_projected_operator(hamiltonian, source, state) == pytest.approx(energies[0], abs=1e-10)
    spin = models.spin()
    raising = model.operator_decomposition(models.located(spin, "Sp", 0)).op
    cross = ed.projected_operator(raising, source, target)
    assert cross.shape == (target.dimension, source.dimension) and cross.nnz > 0
    csr = cross.to_csr()
    assert csr.shape == cross.shape and csr.to_scipy().shape == cross.shape
    np.testing.assert_allclose(cross.to_scipy().toarray(), csr.toarray(), atol=1e-14)
    blocks = ed.symmetry_spectrum_from_local(hamiltonian, spaces.local_spaces, [translation], [0.0], count=2)
    assert len(blocks) == 6 and min(block.eigenvalues[0] for block in blocks) == pytest.approx(energies[0], abs=1e-12)
    with pytest.raises(ValueError):
        ed.prepare_projected(raising, source)


def test_dihedral_family_reproduces_the_sector():
    n = 6
    model = models.heisenberg(lat.chain(n))
    family = [ed.dihedral_irrep(n, 0, reflection=1), ed.dihedral_irrep(n, 0, reflection=-1),
              ed.dihedral_irrep(n, 3, reflection=1), ed.dihedral_irrep(n, 3, reflection=-1),
              ed.dihedral_irrep(n, 1), ed.dihedral_irrep(n, 2)]
    ed.validate_irrep_family(family)
    assert [irrep.dimension for irrep in family] == [1, 1, 1, 1, 2, 2]
    assert family[4].order == 12 and family[4].num_sites == 6 and len(family[4].permutations) == 12
    assert family[4].permutations[0] == list(range(6)) and family[4].matrices[0].shape == (2, 2)
    assert family[4].characters.shape == (12,) and family[4].characters[0] == 2
    with pytest.raises(ValueError):
        ed.validate_irrep_family(family[:3])
    spectrum = ed.finite_spectrum(model, family, sector={"Sz": 0.0})
    assert sum(block.block_dimension for block in spectrum.blocks) == 20
    assert {block.irrep for block in spectrum.blocks} == {irrep.label for irrep in family}
    assert all(block.group_order == 12 for block in spectrum.blocks)
    reference = ed.prepare_product(model, sector={"Sz": 0.0}).eigenstates(count=20, method="dense").eigenvalues
    np.testing.assert_allclose(spectrum.eigenvalues, reference, atol=1e-10)
    problem = ed.prepare_finite(model, family[4], sector={"Sz": 0.0})
    assert problem.basis.finite_group.label == family[4].label and problem.dimension % 2 == 0
    result = problem.eigenstates(count=2)
    hamiltonian = models.compile_hamiltonian(model).op
    np.testing.assert_allclose(result.expectation_values(hamiltonian), result.eigenvalues, atol=1e-10)
    rotation = ed.FiniteGroupIrrep([[1, 2, 3, 4, 5, 0]], [[[np.exp(2j * np.pi / 6)]]], 1, "k1")
    assert rotation.order == 6 and rotation.dimension == 1
    one_momentum = ed.prepare_finite(model, rotation, sector={"Sz": 0.0})
    momentum = ed.prepare_projected(model, sector={"Sz": 0.0}, momentum=[2 * np.pi / 6])
    np.testing.assert_allclose(one_momentum.eigenstates(count=2).eigenvalues, momentum.eigenstates(count=2).eigenvalues, atol=1e-10)
    with pytest.raises(ValueError):
        ed.FiniteGroupIrrep([[1, 0]], [[[2.0]]], 1, "bad")


def test_sparse_assembly_matches_kronecker_products():
    spin = models.spin()
    sz, sp, sm = (spin.op(name).matrix for name in ("Sz", "Sp", "Sm"))
    dimer = models.compile_hamiltonian(
        models.located(spin, "Sz", 0) * models.located(spin, "Sz", 1)
        + 0.5 * (models.located(spin, "Sp", 0) * models.located(spin, "Sm", 1)).add_hc(), [spin, spin])
    states = [[0, 0], [0, 1], [1, 0], [1, 1]]
    expected = np.kron(sz, sz) + 0.5 * (np.kron(sp, sm) + np.kron(sm, sp))
    coo = ed.build_coo_from_local(dimer, states)
    assert coo.shape == (4, 4) and not coo.is_complex
    np.testing.assert_allclose(coo.to_scipy().toarray(), expected, atol=1e-14)
    np.testing.assert_allclose(ed.build_csr_from_local(dimer, states, real=True).toarray(), expected, atol=1e-14)
    np.testing.assert_allclose(coo.to_csr().toarray(), expected, atol=1e-14)
    sector = ed.build_csr_from_local(dimer, [[0, 1], [1, 0]])
    np.testing.assert_allclose(sector.toarray(), expected[1:3, 1:3], atol=1e-14)
    with pytest.raises(ValueError):
        ed.build_coo_from_local(dimer, [[0, 1, 0]])

    basis = ed.enumerate_product_basis([2, 2, 2], [([[0.5, -0.5]] * 3, 0.5)])
    assert basis.dimension == 3 and basis.num_sites == 3
    np.testing.assert_array_equal(basis.states, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    assert len(ed.enumerate_product_basis([2, 3])) == 6

    model = models.heisenberg(lat.chain(8))
    compiled = models.compile_hamiltonian(model)
    sparse = ed.build_sparse_from_local(compiled, constraints=[("Sz", 0.0)])
    assert sparse.dim == 70 and sparse.shape == (70, 70) and sparse.states.shape == (70,)
    sector_energy = ed.lanczos(sparse.to_scipy(), count=1).eigenvalues[0]
    assert sector_energy == pytest.approx(CHAIN_8_GROUND_ENERGY, abs=1e-8)
    assert ed.build_sparse_from_local(compiled).dim == 256

    fermion = models.hubbard(lat.chain(2, bc="open"), t=1.0, U=4.0)
    hamiltonian = models.compile_hamiltonian(fermion).op
    parity = [site.parity for site in models.physical_basis(fermion).sites]
    all_states = [[a, b] for a in range(4) for b in range(4)]
    with_strings = ed.build_csr_from_local(hamiltonian, all_states, parity=parity, real=True).toarray()
    np.testing.assert_allclose(with_strings, ed.prepare_product(fermion).hamiltonian.toarray(), atol=1e-14)
    assert np.linalg.eigvalsh(with_strings)[0] == pytest.approx(ed.prepare_product(fermion).eigenstates(count=1).eigenvalues[0])


def test_thermodynamics_from_a_spectrum():
    energies = np.array([0.0, 1.0, 2.0])
    weights = ed.thermal_weights(energies, 1.0)
    np.testing.assert_allclose(weights, np.exp(-energies) / np.exp(-energies).sum())
    surface = ed.thermal_surface(energies, temperatures=[0.5, 1.0], fields=[0.0])
    assert set(surface) == {"energy", "free_energy", "entropy", "specific_heat", "charge", "susceptibility",
                            "partition", "energy_reference", "log_partition"}
    assert surface["energy"].shape == (2, 1)
    for row, temperature in enumerate([0.5, 1.0]):
        boltzmann = np.exp(-energies / temperature)
        mean = (energies * boltzmann).sum() / boltzmann.sum()
        assert surface["energy"][row, 0] == pytest.approx(mean)
        assert surface["free_energy"][row, 0] == pytest.approx(-temperature * np.log(boltzmann.sum()))
        assert surface["specific_heat"][row, 0] == pytest.approx(
            ((energies ** 2 * boltzmann).sum() / boltzmann.sum() - mean ** 2) / temperature ** 2)
    charged = ed.thermal_surface(energies, temperatures=[1.0], fields=[0.0, 0.5], charges=[1.0, 0.0, -1.0])
    assert charged["charge"].shape == (1, 2) and charged["charge"][0, 1] > charged["charge"][0, 0]
    with pytest.raises(ValueError):
        ed.thermal_weights(energies, 0.0)

    model = models.heisenberg(lat.chain(6))
    spectrum = ed.model_spectrum(model)
    levels = spectrum.eigenvalues
    charges = np.concatenate([[block.quantum_numbers[0]] * len(block) for block in spectrum.blocks])
    order = np.argsort(np.concatenate([block.eigenvalues for block in spectrum.blocks]))
    magnetization = ed.thermal_surface(levels, temperatures=[0.3], fields=[0.0, 1.0], charges=charges[order])
    assert magnetization["charge"][0, 0] == pytest.approx(0.0, abs=1e-12) and magnetization["charge"][0, 1] > 0.5
