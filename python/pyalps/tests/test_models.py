# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""pyalps.models: lattices, local bases, the operator algebra, compiled models,
the Pauli export, the XML bridge and the catalog importer."""

import re

import numpy as np
import pytest

models = pytest.importorskip("pyalps.models")
ed = pytest.importorskip("pyalps.ed")
lat = models.lattice

PAULI = {
    "I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.array([[1, 0], [0, -1]]),
}


def pauli_matrix(pauli):
    """The dense register matrix of a PauliOperator (qubit 0 at the right end)."""
    total = np.zeros((2 ** pauli.num_qubits,) * 2, dtype=complex)
    for string, coefficient in pauli.terms.items():
        matrix = np.array([[1.0]])
        for letter in string:
            matrix = np.kron(matrix, PAULI[letter])
        total += coefficient * matrix
    return total


def ground_energy(model, **sector):
    return ed.prepare_product(model, sector=sector).eigenstates(count=1).eigenvalues[0]


def test_module_surface():
    import pyalps
    assert pyalps.models_c is models._models_c
    from pyalps.models.lattice import chain
    assert chain is lat.chain
    assert "Hamiltonian" in models.__all__ and "chain" in dir(lat)


def test_named_lattices_sites_bonds_and_queries():
    chain = lat.chain(4)
    assert (chain.num_sites, chain.num_bonds, chain.dimension) == (4, 4, 1)
    assert chain.alps_name == "chain lattice" and not chain.is_custom
    assert [(b.source, b.target) for b in chain.bonds] == [(0, 1), (1, 2), (2, 3), (3, 0)]
    assert [b.wraps for b in chain.bonds] == [False, False, False, True]
    assert sorted(chain.neighbors(0)) == [1, 3]
    assert chain.boundary == [lat.Boundary.periodic]
    np.testing.assert_allclose(chain.momenta()[:, 0], [0, np.pi / 2, np.pi, 3 * np.pi / 2])
    assert chain.translation_generators() == [[1, 2, 3, 0]]
    assert chain.periodic_directions() == [0]

    assert lat.chain(4, bc="open").num_bonds == 3
    assert lat.chain(4, bc=lat.Boundary.open).num_bonds == 3
    assert lat.chain(4, bc="open").momenta().shape == (1, 0)
    with pytest.raises(TypeError):
        lat.chain(4, bc="twisted")

    square = lat.square(3, 3)
    assert (square.num_sites, square.num_bonds) == (9, 18)
    assert square.coordinates.shape == (9, 2)
    assert lat.square(3, 3, bc="open").num_bonds == 12
    honeycomb = lat.honeycomb(2, 2)
    assert (honeycomb.num_sites, honeycomb.num_bonds) == (8, 12)
    assert sorted(set(honeycomb.site_types())) == [0, 1]
    assert sorted({b.bond_type for b in honeycomb.bonds}) == [0, 1, 2]
    assert len(honeycomb.sites_by_type(1)) == 4 and len(honeycomb.bonds_by_type(2)) == 4
    assert lat.kagome(2, 2).num_sites == 12
    ladder = lat.ladder(3, bc="open")
    assert (ladder.num_sites, ladder.num_bonds) == (6, 7)
    assert [b.bond_type for b in ladder.bonds].count(1) == 3
    assert lat.simple_cubic(2, 2, 2).num_bonds == 24
    assert lat.triangular(3, 3).num_bonds == 27
    edges = chain.edges()
    assert [(e.source, e.target, e.bond_type) for e in edges] == [(0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 0, 0)]
    assert "chain lattice" in repr(chain)


def test_custom_unit_cells_tile_and_build():
    cell = lat.UnitCell(2)
    assert cell.add_vertex([0.0, 0.0]) == 0
    cell.add_edge(0, 0, 0, [], [1, 0])
    cell.add_edge(0, 0, 0, [], [0, 1])
    assert (cell.dimension, cell.num_vertices, len(cell.edges)) == (2, 1, 2)
    assert cell.edges[1].target_offset == [0, 1]
    periodic = cell.tile([3, 3])
    assert (periodic.num_sites, periodic.num_bonds) == (9, 18) and periodic.is_custom
    cylinder = cell.tile([3, 3], bc=["periodic", "open"])
    assert cylinder.num_bonds == 15
    assert cylinder.boundary == [lat.Boundary.periodic, lat.Boundary.open]
    with pytest.raises(ValueError):
        cell.tile([3, 3], bc=["open"])
    named = lat.Lattice.build(cell, [[1.0, 0.0], [0.0, 1.0]], [2, 2], "open", "square lattice")
    assert named.alps_name == "square lattice" and named.num_bonds == 4
    np.testing.assert_allclose(named.coordinates, [[0, 0], [0, 1], [1, 0], [1, 1]])
    assert named.unit_cell.num_vertices == 1
    assert named.extent == [2, 2] and named.basis_vectors == [[1.0, 0.0], [0.0, 1.0]]


def test_xml_bridge_catalog_and_custom_routes(tmp_path):
    honeycomb = lat.honeycomb(2, 3)
    xml = lat.to_latticegraph_xml(honeycomb, "x")
    assert '<LATTICEGRAPH name="x">' in xml and "<COORDINATE>" in xml
    assert '<EXTENT dimension="2" size="3"/>' in xml
    graph = lat.to_graph_xml(lat.chain(3, bc="open"))
    assert '<GRAPH vertices="3" edges="2">' in graph
    name = lat.custom_lattice_name(honeycomb)
    assert re.fullmatch(r"alps_custom_lattice_[0-9a-f]{16}", name)
    assert name == lat.custom_lattice_name(lat.honeycomb(2, 3))
    assert name != lat.custom_lattice_name(lat.honeycomb(3, 2))
    assert lat.catalog_graph_name(lat.chain(4)) == "chain lattice"
    assert lat.catalog_graph_name(lat.chain(4, bc="open")) == "open chain lattice"
    assert lat.catalog_graph_name(lat.triangular(3, 3, bc="open")) == ""

    square = lat.resolve_lattice(lat.square(2, 3))
    assert square.from_catalog and square.graph == "square lattice"
    assert square.parameters["LATTICE"] == "square lattice"
    assert (square.parameters["L"], square.parameters["W"]) == ("2", "3")
    assert square.parameters["bc1"] == "periodic" and square.library_xml == ""
    assert lat.write_lattice_library(square, str(tmp_path / "unused.xml")) == square.parameters
    assert not (tmp_path / "unused.xml").exists()

    custom = lat.resolve_lattice(honeycomb)
    assert not custom.from_catalog and custom.graph == name
    assert "<UNITCELL" in custom.library_xml
    always = lat.custom_lattice(lat.chain(4))
    assert not always.from_catalog and always.parameters == {"LATTICE": always.graph}
    written = lat.write_lattice_library(custom, str(tmp_path / "honeycomb.xml"))
    assert written["LATTICE"] == name
    assert (tmp_path / "honeycomb.xml").read_text() == custom.library_xml
    assert written["LATTICE_LIBRARY"].endswith("honeycomb.xml")
    with pytest.raises(RuntimeError):
        lat.write_lattice_library(custom, str(tmp_path / "missing" / "dir" / "x.xml"))


def test_built_in_sites_and_custom_sites():
    spin = models.spin()
    assert (spin.dim, spin.statistics) == (2, models.Statistics.bosonic)
    assert {"Sz", "Sp", "Sm", "Sx", "Sy", "S2", "I"} <= set(spin.operator_names)
    np.testing.assert_allclose(spin.op("Sz").matrix, np.diag([0.5, -0.5]))
    np.testing.assert_allclose(spin.op("Sp").matrix, [[0, 1], [0, 0]])
    assert spin.op("Sp").qn_deltas == {"Sz": 1.0}
    assert spin.op("Sp").dag().name == "Sp_dag"
    assert spin.quantum_numbers == {"Sz": [0.5, -0.5]} and spin.parity_operator == [1.0, 1.0]
    assert models.spin(1.0).dim == 3 and models.boson(2).dim == 3
    fermion = models.spinful_fermion()
    assert (fermion.dim, fermion.statistics) == (4, models.Statistics.fermionic)
    assert fermion.parity_operator == [1.0, -1.0, -1.0, 1.0]
    assert fermion.op("c_up").parity == models.Parity.odd and fermion.op("n_up").parity == models.Parity.even
    assert models.spinless_fermion().op("c")(0, 1) == 1.0

    custom = models.SiteBasis(
        2, {"X": [[0, 1], [1, 0]], "Z": [[1, 0], [0, -1]], "N": np.diag([0.0, 1.0])},
        quantum_numbers={"N": [0.0, 1.0]}, labels=["empty", "full"])
    assert custom.has_op("I") and custom.op("I").is_identity()
    assert custom.op("Z").qn_deltas == {"N": 0.0} and custom.op("X").qn_deltas == {}
    assert custom.basis_labels == ["empty", "full"]
    assert models.compute_qn_deltas([[0, 0], [1, 0]], {"N": [0.0, 1.0]}) == {"N": 1.0}
    odd = models.SiteBasis(2, {"c": [[0, 1], [0, 0]]}, quantum_numbers={"N": [0.0, 1.0]},
                           statistics="fermionic", parity_operator=[1.0, -1.0], odd_operators=["c"])
    assert odd.op("c").parity == models.Parity.odd
    with pytest.raises(ValueError):
        models.SiteBasis(2, {"X": np.eye(3)})



def test_catalog_spellings_and_named_adjoints():
    spin, fermion = models.spin(), models.spinful_fermion()
    assert models.catalog_operator_names["Sp"] == "Splus" and models.catalog_operator_names["n_dn"] == "n_down"
    assert spin.op("Splus").name == "Sp" and spin.has_op("Sminus") and not spin.has_op("cdag_up")
    assert fermion.op("cdag_up").name == "c_dag_up" and fermion.op("double_occupancy").name == "n_updn"
    assert models.located(fermion, "c_down", 1).terms[0].factors[0].op.name == "c_dn"
    assert spin.adjoint(spin.op("Sp")).name == "Sm" and spin.adjoint(spin.op("Sz")).name == "Sz"
    assert fermion.adjoint(fermion.op("c_up")).name == "c_dag_up"
    hop = (models.located(fermion, "c_dag_up", 0) * models.located(fermion, "c_up", 2)).dag()
    assert [(f.op.name, f.site) for f in hop.terms[0].factors] == [("c_dag_up", 2), ("c_up", 0)]
    flip = (models.located(spin, "Sp", 0) * models.located(spin, "Sm", 1)).add_hc()
    assert [f.op.name for f in flip.terms[1].factors] == ["Sp", "Sm"]
    with pytest.raises(IndexError):
        spin.op("cdag_up")

    # The other direction: kit names on a catalog model.
    chain = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 4, "J": 1.0})
    np.testing.assert_array_equal(chain.compile_operator("Sp", site=2).terms[0].factors[0].matrix,
                                  chain.compile_operator("Splus", site=2).terms[0].factors[0].matrix)
    hubbard = models.CatalogModel({"MODEL": "fermion Hubbard", "LATTICE": "open chain lattice", "L": 3,
                                   "t": 1.0, "U": 4.0})
    assert len(hubbard.compile_operator("n_updn")) == len(hubbard.compile_operator("double_occupancy")) == 3


def test_site_operator_sugar_observables_and_matrices():
    spin, fermion = models.spin(), models.spinful_fermion()
    assert isinstance(spin.Sz, models.LocatedOperator) and spin.Sz.name == "Sz"
    for built in (spin.Sp(2), spin["Sp"](2), spin.Splus(2)):
        factor = built.terms[0].factors[0]
        assert (factor.op.name, factor.site) == ("Sp", 2)
    assert fermion.c_dag_up(1).terms[0].factors[0].op.name == "c_dag_up"
    assert spin.dim == 2  # attributes still win over operator names
    with pytest.raises(AttributeError):
        spin.cdag_up
    with pytest.raises(KeyError):
        spin["no_such_operator"]

    model = models.heisenberg(lat.chain(6))
    compiled = models.compile_hamiltonian(model).op
    assert models.compile_observable(model, compiled) is compiled
    total = models.compile_observable(model, "Sz")
    assert len(total) == 6 and [term.factors[0].site for term in total.terms] == list(range(6))
    pair = models.compile_observable(model, spin.Sz(0) * spin.Sz(3))
    assert len(pair) == 1 and pair.site_dims == [2] * 6
    assert len(models.compile_observable(models.physical_basis(model), spin.Sz(1))) == 1
    with pytest.raises(ValueError):
        models.compile_observable(model, "c_up")
    with pytest.raises(ValueError):
        models.compile_observable(model, models.spin(1.0).Sz(0))  # another local space
    with pytest.raises(TypeError):
        models.compile_observable(models.physical_basis(model), "Sz")

    dense = model.to_dense_matrix(sector={"Sz": 0.0})
    assert dense.shape == (20, 20)
    np.testing.assert_allclose(np.linalg.eigvalsh(dense)[0], ground_energy(model, Sz=0.0), atol=1e-12)
    assert model.to_sparse_matrix().shape == (64, 64)
    catalog = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 6, "J": 1.0,
                                   "CONSERVED_QUANTUMNUMBERS": "Sz", "Sz_total": 0})
    np.testing.assert_allclose(np.linalg.eigvalsh(catalog.to_dense_matrix()), np.linalg.eigvalsh(dense), atol=1e-12)
    assert len(models.compile_observable(catalog, "Sz")) == 6

def test_term_algebra_and_compilation():
    spin = models.spin()
    sites = [spin, spin]
    hopping = models.located(spin, "Sp", 0) * models.located(spin, "Sm", 1)
    exchange = 0.5 * hopping.add_hc() + models.located(spin, "Sz", 0) * models.located(spin, "Sz", 1)
    assert len(exchange) == 3 and not exchange.empty
    assert (exchange - exchange).canonicalise().empty
    doubled = (2 * exchange).canonicalise()
    assert [t.effective_coeff for t in doubled.terms] == [1.0, 1.0, 2.0]
    assert exchange.dag().canonicalise().terms[0].effective_coeff == 0.5
    assert (exchange / 2).canonicalise().terms[2].effective_coeff == 0.5
    assert [f.op.name for f in hopping.terms[0].factors] == ["Sp", "Sm"]
    assert hopping.terms[0].qn_delta() == {"Sz": 0.0}

    compiled = models.compile_hamiltonian(exchange, sites, conserved_quantum_numbers=["Sz"])
    assert compiled.site_dims == [2, 2] and len(compiled) == 3
    assert compiled.terms[0].sites == [0, 1] and compiled.terms[0].jordan_wigner == []
    np.testing.assert_allclose(compiled.terms[2].factors[0].matrix, np.diag([0.5, -0.5]))
    observable = models.compile_operator(models.located(spin, "Sp", 0), sites)
    assert len(observable) == 1 and observable.terms[0].factors[0].site == 0
    with pytest.raises(ValueError):
        models.validate_hermiticity(models.located(spin, "Sp", 0), sites)
    with pytest.raises(ValueError):
        models.compile_hamiltonian(models.located(spin, "Sp", 0), sites, conserved_quantum_numbers=["Sz"])

    # Compiled operators can also be built directly.
    sz = spin.op("Sz").matrix
    term = models.LocalTerm(1.0, [models.LocalFactor(0, sz), models.LocalFactor(1, sz)])
    op = models.CompiledOperator([2, 2], [term])
    assert op.num_sites == 2 and op.terms[0].factors[1].dim == 2
    assert op == models.CompiledOperator([2, 2], [term]) and op != compiled


def test_in_place_sums_append_in_order_and_change_the_object():
    spin = models.spin()
    a = 0.5 * models.located(spin, "Sp", 0) * models.located(spin, "Sm", 1)
    b = models.located(spin, "Sz", 0) - 2j * models.located(spin, "Sz", 1)
    reference = a + b
    alias = total = models.TermSum() + a
    total += b
    assert total is alias and len(total) == len(reference) == 3
    assert [t.effective_coeff for t in total.terms] == [t.effective_coeff for t in reference.terms]
    assert [[f.op.name for f in t.factors] for t in total.terms] == [["Sp", "Sm"], ["Sz"], ["Sz"]]
    assert len(a) == 1  # the right operand is untouched
    total -= b
    assert (total - a).canonicalise().empty
    total += total
    assert len(total) == 10


def test_hamiltonian_templates_reproduce_the_built_in_model():
    spin = models.spin()
    chain = lat.chain(4)
    model = models.Hamiltonian(chain, spin, conserved_quantum_numbers=["Sz"])
    sz_sz = lambda i, j: models.located(spin, "Sz", i) * models.located(spin, "Sz", j)  # noqa: E731
    flip = lambda i, j: models.located(spin, "Sp", i) * models.located(spin, "Sm", j)  # noqa: E731
    assert model.add_bond(1.0, sz_sz).add_bond(0.5, flip, add_hc=True) is model
    assert model.num_sites == 4 and model.site_dims == [2, 2, 2, 2]
    assert model.conserved_quantum_numbers == ["Sz"]
    assert model.fermion_convention == models.FermionConvention.strict
    assert len(model.expand()) == 12
    compiled = model.compile()
    assert len(compiled.op) == 12 and len(compiled.edges) == 4 and compiled.site_types == [0] * 4
    assert compiled.diag_observables == {0: {"Sz": [0.5, -0.5]}} and compiled.spin_basis == {}
    compiled.spin_basis = {0: [0.5, -0.5]}  # sampler metadata is writable
    assert compiled.spin_basis == {0: [0.5, -0.5]}
    assert ground_energy(model, Sz=0.0) == pytest.approx(-2.0, abs=1e-12)
    assert ground_energy(models.heisenberg(chain), Sz=0.0) == pytest.approx(-2.0, abs=1e-12)
    assert model.site_for(2).dim == 2 and model.lattice.num_bonds == 4
    field = models.Hamiltonian(chain, spin).add_site(-1.0, lambda i: models.located(spin, "Sz", i))
    assert ground_energy(field) == pytest.approx(-2.0, abs=1e-12)
    assert field.physical_basis().conserved_quantum_numbers == []

    mixed = models.Hamiltonian(lat.honeycomb(1, 2), {0: models.spin(0.5), 1: models.spin(1.0)})
    assert mixed.site_dims == [2, 3, 2, 3]
    hubbard = models.hubbard(lat.chain(2, bc="open"), t=1.0, U=4.0)
    assert hubbard.site_dims == [4, 4] and hubbard.conserved_quantum_numbers == ["N", "Sz"]
    basis = models.physical_basis(hubbard)
    assert basis.site_dims == [4, 4] and basis.sites[0].statistics == models.Statistics.fermionic
    assert basis.sites[1].quantum_numbers["N"] == [0.0, 1.0, 1.0, 2.0]
    assert basis == hubbard.physical_basis() and basis != models.physical_basis(model)
    for factory in (models.ising, models.xxz, models.xy, models.transverse_field_ising, models.heisenberg_d,
                    models.biquadratic, models.j1j2, models.t_j, models.hardcore_boson_model, models.bose_hubbard,
                    models.spinless_fermion_model):
        assert factory(chain).num_sites == 4
    assert models.kondo_lattice(lat.honeycomb(1, 2)).site_dims == [4, 2, 4, 2]  # two vertex types
    assert models.j1j2(lat.nnn_chain(4)).lattice.num_bonds == 8
    observable = model.operator_decomposition(models.located(spin, "Sp", 1))
    assert len(observable.op) == 1 and observable.edges == compiled.edges


def test_python_specifications_are_evaluated_when_added():
    import gc
    import weakref

    spin = models.spin()
    model = models.Hamiltonian(lat.chain(4), spin)
    spec = lambda i, j: models.located(spin, "Sz", i) * models.located(spin, "Sz", j)  # noqa: E731
    alive = weakref.ref(spec)
    model.add_bond(1.0, spec)
    del spec
    gc.collect()
    assert alive() is None  # the model keeps no reference to the callable
    assert len(model.expand()) == 4
    calls = []
    model.add_site(-0.5, lambda i: calls.append(i) or models.located(spin, "Sz", i))
    assert calls == [0, 1, 2, 3] and len(model.expand()) == 8
    model.compile()
    assert calls == [0, 1, 2, 3]  # not called again
    assert model.add_bond(1.0, spec_never_called := (lambda i, j: 1 / 0), bond_type=7) is model  # no such bonds
    del spec_never_called
    with pytest.raises(ZeroDivisionError):
        model.add_site(1.0, lambda i: 1 / 0)
    with pytest.raises(TypeError, match="must return a TermSum, not str"):
        model.add_site(1.0, lambda i: "not a TermSum")


def test_physical_basis_values():
    basis = models.physical_basis([models.spin(), models.spinless_fermion()], conserved_quantum_numbers=["N"])
    assert basis.num_sites == 2 and basis.conserved_quantum_numbers == ["N"]
    assert basis.sites[0].labels == ["0", "1"] and basis.sites[0].quantum_numbers["N"] == [0.0, 0.0]
    assert basis.sites[1].parity == [1.0, -1.0] and basis.sites[1].statistics == models.Statistics.fermionic
    models.validate_basis(basis)
    rebuilt = models.PhysicalBasis(basis.sites, ["N"], "strict")
    assert rebuilt == basis
    local = models.LocalBasis(["up", "down"], "bosonic", {"Sz": [0.5, -0.5]}, [1.0, 1.0])
    assert local.dim == 2 and local != basis.sites[0]
    with pytest.raises(ValueError):
        models.validate_basis(models.PhysicalBasis([models.LocalBasis(["a"], "bosonic", {}, [2.0])]))
    with pytest.raises(ValueError):
        models.physical_basis([models.spin()], conserved_quantum_numbers=["N"])


def test_pauli_export_and_qubit_encoding():
    model = models.heisenberg(lat.chain(4, bc="open"))
    pauli = models.export_pauli(model)
    assert pauli.num_qubits == 4 and pauli.site_qubits == [1, 1, 1, 1]
    assert len(pauli) == 9 and all(abs(c - 0.25) < 1e-15 for c in pauli.terms.values())
    assert {s[-2:] for s in pauli.terms if s.startswith("II")} == {"XX", "YY", "ZZ"}
    assert pauli.state_index([0, 1, 0, 1]) == 10
    assert pauli.discarded_l1 == 0 and pauli.fermion_mode_qubits == []
    assert pauli.coefficient_l1() == pytest.approx(2.25)
    register = pauli_matrix(pauli)
    reference = ed.prepare_product(model).hamiltonian.toarray()
    np.testing.assert_allclose(np.linalg.eigvalsh(register), np.linalg.eigvalsh(reference), atol=1e-12)

    basis = models.physical_basis(model)
    gray = models.QubitEncoding(basis, scheme="gray")
    assert gray.scheme == models.QubitEncodingScheme.gray and gray.num_qubits == 4
    assert gray.site_offsets == [0, 1, 2, 3]
    assert gray.decode_basis_bits(gray.encode_basis_bits([1, 0, 1, 0])) == [1, 0, 1, 0]
    assert models.make_qubit_encoding(model).state_index([1, 1, 1, 1]) == 15
    physical = np.zeros(16)
    physical[5] = 1.0
    encoded = gray.encode_statevector(physical)
    decoded, discarded = gray.decode_statevector(encoded)
    np.testing.assert_array_equal(decoded, physical)
    assert discarded == 0 and gray.physical_weight(encoded) == pytest.approx(1.0)
    same = models.export_pauli(models.compile_hamiltonian(model).op, basis, encoding=gray)
    assert same.terms == pauli.terms
    with pytest.raises(ValueError):
        models.export_pauli(models.compile_hamiltonian(model).op, models.physical_basis([models.spin()] * 3))

    spin = models.spin()
    raising = model.operator_decomposition(models.located(spin, "Sp", 2)).op
    with pytest.raises(ValueError):
        models.export_pauli(raising, basis)
    exported = models.export_pauli(raising, basis, require_hermitian=False, require_even=False,
                                   require_conservation=False)
    assert set(exported.terms) == {"IXII", "IYII"}

    hubbard = models.hubbard(lat.chain(2, bc="open"), t=1.0, U=4.0)
    fermions = models.export_pauli(hubbard)
    assert fermions.num_qubits == 4 and fermions.fermion_mode_qubits == [1, 0, 3, 2]
    hubbard_reference = ed.prepare_product(hubbard).hamiltonian.toarray()
    np.testing.assert_allclose(np.linalg.eigvalsh(pauli_matrix(fermions)), np.linalg.eigvalsh(hubbard_reference), atol=1e-12)

    boson = models.bose_hubbard(lat.chain(2, bc="open"), n_max=2)  # dimension 3 per site: one unused code word
    encoding = models.make_qubit_encoding(boson)
    assert encoding.site_qubits == [2, 2]
    projector = models.local_valid_projector(encoding, 0)
    penalty = models.invalid_state_penalty(encoding, 10.0)
    assert len(projector) > 0 and len(penalty) > 0
    assert np.allclose(np.linalg.eigvalsh(pauli_matrix(projector)), [0] * 4 + [1] * 12)


def test_catalog_model_matches_the_kit():
    chain = lat.chain(4)
    catalog = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 4, "J": 1.0,
                                   "local_S": 0.5, "CONSERVED_QUANTUMNUMBERS": "Sz"})
    assert catalog.num_sites == 4 and catalog.site_dims == [2, 2, 2, 2]
    assert catalog.conserved_quantum_numbers() == ["Sz"] and catalog.sector() == {}
    assert catalog.local_parity() == [[1.0, 1.0]] * 4
    assert catalog.local_charge("Sz") == [[-0.5, 0.5]] * 4  # the XML spin basis counts Sz upwards
    assert [(e.source, e.target) for e in catalog.edges] == [(0, 1), (1, 2), (2, 3), (3, 0)]
    assert catalog.parameters["L"] == "4"
    hamiltonian = catalog.compile_hamiltonian()
    assert hamiltonian.site_dims == [2, 2, 2, 2]
    basis = catalog.physical_basis()
    assert basis.conserved_quantum_numbers == ["Sz"] and basis.sites[0].labels == ["S=0.5 Sz=-0.5", "S=0.5 Sz=0.5"]
    spaces = ed.projected_spaces(basis)
    problem = ed.prepare_projected(hamiltonian, spaces, sector={"Sz": 0.0})
    assert problem.dimension == 6
    assert problem.eigenstates(count=1).eigenvalues[0] == pytest.approx(-2.0, abs=1e-12)
    assert ground_energy(models.heisenberg(chain), Sz=0.0) == pytest.approx(-2.0, abs=1e-12)
    sz1 = catalog.compile_operator("Sz", site=1)
    assert len(sz1) == 1 and sz1.terms[0].factors[0].site == 1
    assert len(catalog.compile_operator("Sz")) == 4
    with pytest.raises(ValueError):
        catalog.compile_operator("Sz", site=1, bond=0)
    with pytest.raises(IndexError):
        catalog.compile_operator("Sz", site=9)
    catalog.set_parameters({"Sz_total": 1})
    assert catalog.sector() == {"Sz": 1.0}
    assert "spin" in repr(catalog)

    # The same catalog model on a kit lattice, against the kit Hamiltonian with per-type couplings.
    ladder = lat.ladder(3, bc="open")
    on_ladder = models.CatalogModel(ladder, {"MODEL": "spin", "J0": 1.0, "J1": 0.4, "local_S": 0.5,
                                             "CONSERVED_QUANTUMNUMBERS": "Sz"})
    assert on_ladder.num_sites == 6 and on_ladder.num_bonds == 7
    spin = models.spin()
    kit = models.heisenberg(ladder, J=1.0, bond_type=0)
    kit.add_bond(0.4, lambda i, j: (models.located(spin, "Sz", i) * models.located(spin, "Sz", j)
                                    + 0.5 * (models.located(spin, "Sp", i) * models.located(spin, "Sm", j)).add_hc()),
                 bond_type=1)
    catalog_problem = ed.prepare_projected(on_ladder.compile_hamiltonian(), ed.projected_spaces(on_ladder.physical_basis()),
                                           sector={"Sz": 0.0})
    assert catalog_problem.eigenstates(count=1).eigenvalues[0] == pytest.approx(ground_energy(kit, Sz=0.0), abs=1e-11)
    with pytest.raises(RuntimeError):
        models.CatalogModel({"MODEL": "spin", "LATTICE": "no such lattice", "L": 4})
    with pytest.raises(TypeError):
        models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": [1, 2]})
