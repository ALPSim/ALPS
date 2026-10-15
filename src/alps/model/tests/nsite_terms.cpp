/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// n-site terms: <PLAQUETTE> in lattices, <PLAQUETTETERM> in models, and their
// contribution to hamiltonian_matrix.
//
// Reference energies: the ladder values come from an independent Python
// exact-diagonalization prototype (dense Jordan-Wigner, mode order 2*site+spin),
// the "mode U" spectrum is analytic.

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

const char* fermion_library = R"(
<SITEBASIS name="fermion">
  <QUANTUMNUMBER name="Nup" min="0" max="1" type="fermionic"/>
  <QUANTUMNUMBER name="Ndown" min="0" max="1" type="fermionic"/>
  <OPERATOR name="c_down" matrixelement="1"><CHANGE quantumnumber="Ndown" change="-1"/></OPERATOR>
  <OPERATOR name="cdag_down" matrixelement="1"><CHANGE quantumnumber="Ndown" change="1"/></OPERATOR>
  <OPERATOR name="c_up" matrixelement="1"><CHANGE quantumnumber="Nup" change="-1"/></OPERATOR>
  <OPERATOR name="cdag_up" matrixelement="1"><CHANGE quantumnumber="Nup" change="1"/></OPERATOR>
  <OPERATOR name="n" matrixelement="Nup+Ndown"/>
  <OPERATOR name="n_up" matrixelement="Nup"/>
  <OPERATOR name="n_down" matrixelement="Ndown"/>
</SITEBASIS>
<BASIS name="fermion">
  <SITEBASIS ref="fermion"/>
  <CONSTRAINT quantumnumber="Nup" value="Nup_total"/>
  <CONSTRAINT quantumnumber="Ndown" value="Ndown_total"/>
</BASIS>
<BONDOPERATOR name="fermion_hop" source="x" target="y">
  cdag_up(x)*c_up(y)+cdag_up(y)*c_up(x)+cdag_down(x)*c_down(y)+cdag_down(y)*c_down(x)
</BONDOPERATOR>
)";

const char* ladder_cell = R"(
<LATTICE name="chain lattice" dimension="1">
  <BASIS><VECTOR>1</VECTOR></BASIS>
</LATTICE>
<UNITCELL name="ladder cell" dimension="1" vertices="2">
  <VERTEX id="1"/><VERTEX id="2"/>
  <EDGE type="1"><SOURCE vertex="1" offset="0"/><TARGET vertex="2" offset="0"/></EDGE>
  <EDGE type="0"><SOURCE vertex="1" offset="0"/><TARGET vertex="1" offset="1"/></EDGE>
  <EDGE type="0"><SOURCE vertex="2" offset="0"/><TARGET vertex="2" offset="1"/></EDGE>
  <PLAQUETTE type="0">
    <SITE vertex="1" offset="0"/><SITE vertex="2" offset="0"/>
    <SITE vertex="2" offset="1"/><SITE vertex="1" offset="1"/>
  </PLAQUETTE>
</UNITCELL>
<LATTICEGRAPH name="open ladder">
  <FINITELATTICE><LATTICE ref="chain lattice"/><EXTENT dimension="1" size="L"/>
    <BOUNDARY type="open"/></FINITELATTICE>
  <UNITCELL ref="ladder cell"/>
</LATTICEGRAPH>
<LATTICEGRAPH name="periodic ladder">
  <FINITELATTICE><LATTICE ref="chain lattice"/><EXTENT dimension="1" size="L"/>
    <BOUNDARY type="periodic"/></FINITELATTICE>
  <UNITCELL ref="ladder cell"/>
</LATTICEGRAPH>
)";

const char* hop = "-t*fermion_hop(i,j)";
const char* pair_hop =
  "-K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))"
  "   + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )";

struct Files {
  std::filesystem::path dir;
  Files(const std::string& lattices, const std::string& hamiltonians) {
    static int counter = 0;
    dir = std::filesystem::temp_directory_path()
      / ("alps_nsite_" + std::to_string(::getpid()) + "_" + std::to_string(counter++));
    std::filesystem::create_directories(dir);
    std::ofstream(dir / "lattices.xml") << "<LATTICES>" << lattices << "</LATTICES>\n";
    std::ofstream(dir / "models.xml") << "<MODELS>" << fermion_library << hamiltonians << "</MODELS>\n";
  }
  ~Files() { std::error_code ec; std::filesystem::remove_all(dir, ec); }
  alps::Parameters parms() const {
    alps::Parameters p;
    p["LATTICE_LIBRARY"] = (dir / "lattices.xml").string();
    p["MODEL_LIBRARY"] = (dir / "models.xml").string();
    return p;
  }
};

std::string hubbard_model(const std::string& name, const std::string& plaquette_term) {
  return "<HAMILTONIAN name=\"" + name + "\">"
         "<PARAMETER name=\"t\" default=\"1\"/><PARAMETER name=\"U\" default=\"0\"/>"
         "<PARAMETER name=\"K\" default=\"0\"/><BASIS ref=\"fermion\"/>"
         "<SITETERM site=\"i\">U*n_up(i)*n_down(i)</SITETERM>"
         "<BONDTERM source=\"i\" target=\"j\">" + std::string(hop) + "</BONDTERM>"
         + plaquette_term + "</HAMILTONIAN>";
}

// Eigenvalues of a small real symmetric matrix (cyclic Jacobi).
std::vector<double> eigenvalues(std::vector<std::vector<double> > a) {
  const std::size_t n = a.size();
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q)
        off += a[p][q] * a[p][q];
    if (off < 1e-26)
      break;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q) {
        if (std::abs(a[p][q]) < 1e-300)
          continue;
        double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        double t = (theta >= 0 ? 1 : -1) / (std::abs(theta) + std::sqrt(theta * theta + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
  }
  std::vector<double> e(n);
  for (std::size_t i = 0; i < n; ++i)
    e[i] = a[i][i];
  std::sort(e.begin(), e.end());
  return e;
}

std::vector<double> spectrum(const alps::Parameters& p) {
  alps::hamiltonian_matrix<alps::numeric::matrix<double> > h(p);
  const alps::numeric::matrix<double>& m = h.matrix();
  std::vector<std::vector<double> > a(h.dimension(), std::vector<double>(h.dimension()));
  for (std::size_t i = 0; i < h.dimension(); ++i)
    for (std::size_t j = 0; j < h.dimension(); ++j)
      a[i][j] = m(i, j);
  for (std::size_t i = 0; i < a.size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
      EXPECT_NEAR(a[i][j], a[j][i], 1e-12) << "Hamiltonian not symmetric";
  return eigenvalues(a);
}

alps::Parameters sector(alps::Parameters p, int nup, int ndown) {
  p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown";
  p["Nup_total"] = nup;
  p["Ndown_total"] = ndown;
  return p;
}

} // namespace

TEST(NSiteLattice, UnitCellPlaquettesFollowBoundaries) {
  Files f(ladder_cell, "");
  alps::Parameters p = f.parms();
  p["L"] = 3;
  p["LATTICE"] = "open ladder";
  alps::graph_helper<> open(p);
  ASSERT_EQ(open.num_plaquettes(), 2u);   // the plaquette leaving the lattice is dropped
  EXPECT_EQ(open.plaquettes()[0].sites, (std::vector<std::size_t>{0, 1, 3, 2}));
  EXPECT_EQ(open.plaquettes()[1].sites, (std::vector<std::size_t>{2, 3, 5, 4}));
  p["LATTICE"] = "periodic ladder";
  alps::graph_helper<> periodic(p);
  ASSERT_EQ(periodic.num_plaquettes(), 3u);
  EXPECT_EQ(periodic.plaquettes()[2].sites, (std::vector<std::size_t>{4, 5, 1, 0}));
  EXPECT_EQ(periodic.num_bonds(), 9u);    // the bonds are unchanged
}

TEST(NSiteLattice, PeriodicPlaquetteOnTooSmallLatticeThrows) {
  Files f(ladder_cell, "");
  alps::Parameters p = f.parms();
  p["L"] = 1;
  p["LATTICE"] = "periodic ladder";
  EXPECT_THROW(alps::graph_helper<> g(p), std::runtime_error);
}

TEST(NSiteLattice, GraphPlaquettes) {
  Files f("<GRAPH name=\"g\" vertices=\"3\"><EDGE source=\"1\" target=\"2\"/>"
          "<PLAQUETTE type=\"2\"><SITE id=\"3\"/><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE></GRAPH>", "");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "g";
  alps::graph_helper<> g(p);
  ASSERT_EQ(g.num_plaquettes(), 1u);
  EXPECT_EQ(g.plaquettes()[0].type, 2u);
  EXPECT_EQ(g.plaquettes()[0].sites, (std::vector<std::size_t>{2, 0, 1}));
}

TEST(NSiteLattice, GraphPlaquettesSurviveLibraryRoundTrip) {
  Files f("<GRAPH name=\"g\" vertices=\"3\"><EDGE source=\"1\" target=\"2\"/>"
          "<PLAQUETTE type=\"2\"><SITE id=\"3\"/><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE>"
          "<PLAQUETTE><SITE id=\"2\"/><SITE id=\"3\"/></PLAQUETTE></GRAPH>", "");
  alps::LatticeLibrary lib(f.parms());
  std::ostringstream out;
  {
    alps::oxstream xml(out);
    lib.write_xml(xml);
  }
  std::istringstream in(out.str());
  alps::LatticeLibrary reread;
  reread.read_xml(in);
  const alps::plaquette_list& p = reread.graph_plaquettes("g");
  ASSERT_EQ(p.size(), 2u);
  EXPECT_EQ(p[0].type, 2u);
  EXPECT_EQ(p[0].sites, (std::vector<std::size_t>{2, 0, 1}));
  EXPECT_EQ(p[1].type, 0u);
  EXPECT_EQ(p[1].sites, (std::vector<std::size_t>{1, 2}));
}

TEST(NSiteLattice, GraphPlaquetteWithRepeatedSiteThrows) {
  Files f("<GRAPH name=\"bad\" vertices=\"2\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"1\"/></PLAQUETTE></GRAPH>", "");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "bad";
  EXPECT_THROW(alps::graph_helper<> g(p), std::runtime_error);
}

TEST(NSiteHamiltonian, LadderPairHoppingMatchesReference) {
  struct Case { const char* lattice; int L; double K; int nup, ndown; double e0; };
  const Case cases[] = {
    {"open ladder", 3, 0.7, 2, 2, -5.3669606044},
    {"open ladder", 3, -1.3, 2, 2, -5.0797988157},
    {"periodic ladder", 3, 0.9, 3, 2, -5.0117188128},
  };
  Files f(ladder_cell, hubbard_model("pair hopping",
          "<PLAQUETTETERM sites=\"i j k l\">" + std::string(pair_hop) + "</PLAQUETTETERM>"));
  for (const Case& c : cases) {
    alps::Parameters p = sector(f.parms(), c.nup, c.ndown);
    p["LATTICE"] = c.lattice; p["L"] = c.L; p["MODEL"] = "pair hopping";
    p["U"] = 8.0; p["K"] = c.K;
    EXPECT_NEAR(spectrum(p)[0], c.e0, 1e-9) << c.lattice << " K=" << c.K;
  }
}

TEST(NSiteHamiltonian, ZeroCouplingReproducesSiteAndBondModel) {
  Files f(ladder_cell, hubbard_model("with K", "<PLAQUETTETERM sites=\"i j k l\">" + std::string(pair_hop) + "</PLAQUETTETERM>")
                       + hubbard_model("without K", ""));
  alps::Parameters p = sector(f.parms(), 2, 1);
  p["LATTICE"] = "open ladder"; p["L"] = 3; p["U"] = 4.0; p["K"] = 0.0;
  p["MODEL"] = "with K";
  std::vector<double> a = spectrum(p);
  p["MODEL"] = "without K";
  std::vector<double> b = spectrum(p);
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i)
    EXPECT_NEAR(a[i], b[i], 1e-12);
}

// Open 4-site chain plus U n_{a,up} n_{a,down} on the hopping eigenmode
// a = sum_i phi_i c_i. In the site basis this is one 4-site term with pieces on
// 1, 2, 3 and 4 distinct sites; the spectrum is
//   sum of occupied mode energies + U [mode a occupied by both spins].
// The plaquette lists the sites in a scrambled order to exercise the signs.
TEST(NSiteHamiltonian, InteractionOnDelocalizedModeIsAnalytic) {
  const int L = 4, ka = 1;
  const double U = 3.7, pi = std::acos(-1.0);
  const std::vector<int> order = {2, 0, 3, 1};
  std::vector<double> eps, phi;
  for (int k = 0; k < L; ++k)
    eps.push_back(-2 * std::cos(pi * (k + 1) / (L + 1)));
  for (int i = 0; i < L; ++i)
    phi.push_back(std::sqrt(2.0 / (L + 1)) * std::sin(pi * (ka + 1) * (i + 1) / (L + 1)));
  auto lin = [&](const std::string& op) {
    std::ostringstream s;
    s.precision(17);
    s << "(";
    for (int r = 0; r < L; ++r)
      s << (r ? " + " : "") << "(" << phi[order[r]] << ")*" << op << "(" << "ijkl"[r] << ")";
    return s.str() + ")";
  };
  std::string term = "U*" + lin("cdag_up") + "*" + lin("c_up") + "*" + lin("cdag_down") + "*" + lin("c_down");
  std::ostringstream graph;
  graph << "<GRAPH name=\"chain4\" vertices=\"4\">";
  for (int i = 1; i < L; ++i)
    graph << "<EDGE source=\"" << i << "\" target=\"" << i + 1 << "\"/>";
  graph << "<PLAQUETTE>";
  for (int s : order)
    graph << "<SITE id=\"" << s + 1 << "\"/>";
  graph << "</PLAQUETTE></GRAPH>";
  Files f(graph.str(), "<HAMILTONIAN name=\"mode U\"><PARAMETER name=\"t\" default=\"1\"/>"
          "<PARAMETER name=\"U\" default=\"0\"/><BASIS ref=\"fermion\"/>"
          "<BONDTERM source=\"i\" target=\"j\">" + std::string(hop) + "</BONDTERM>"
          "<PLAQUETTETERM sites=\"i j k l\">" + term + "</PLAQUETTETERM></HAMILTONIAN>");
  for (auto [nup, ndown] : std::vector<std::pair<int, int> >{{2, 1}, {3, 2}}) {
    std::vector<double> expected;
    for (int up = 0; up < (1 << L); ++up)
      for (int dn = 0; dn < (1 << L); ++dn) {
        if (__builtin_popcount(up) != nup || __builtin_popcount(dn) != ndown)
          continue;
        double e = 0;
        for (int k = 0; k < L; ++k)
          e += eps[k] * (((up >> k) & 1) + ((dn >> k) & 1));
        expected.push_back(e + U * (((up >> ka) & 1) && ((dn >> ka) & 1)));
      }
    std::sort(expected.begin(), expected.end());
    alps::Parameters p = sector(f.parms(), nup, ndown);
    p["GRAPH"] = "chain4"; p["MODEL"] = "mode U"; p["U"] = U;
    std::vector<double> e = spectrum(p);
    ASSERT_EQ(e.size(), expected.size());
    for (std::size_t i = 0; i < e.size(); ++i)
      EXPECT_NEAR(e[i], expected[i], 1e-10) << "sector " << nup << "," << ndown << " level " << i;
  }
}

TEST(NSiteHamiltonian, TermLeavingTheSectorThrows) {
  Files f("<GRAPH name=\"g\" vertices=\"4\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/><SITE id=\"3\"/><SITE id=\"4\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"pairing\"><BASIS ref=\"fermion\"/><PLAQUETTETERM sites=\"i j k l\">"
          "cdag_up(i)*cdag_down(j)*cdag_up(k)*cdag_down(l) + c_down(l)*c_up(k)*c_down(j)*c_up(i)"
          "</PLAQUETTETERM></HAMILTONIAN>");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "g"; p["MODEL"] = "pairing";
  p["CONSERVED_QUANTUMNUMBERS"] = "";
  EXPECT_EQ(spectrum(p).size(), 256u);
  alps::Parameters q = sector(p, 1, 1);
  EXPECT_THROW(spectrum(q), std::runtime_error);
}

// Site, bond and plaquette terms that change a quantum number fixed by the
// sector stop the run, with or without momentum sectors. Elements leaving the
// sector used to be dropped silently, which removed the pairing from the
// spectrum.
std::string build_error(const alps::Parameters& p) {
  try {
    alps::hamiltonian_matrix<alps::numeric::matrix<std::complex<double> > > h(p);
    h.matrix();
  } catch (std::runtime_error& e) {
    return e.what();
  }
  return "";
}

const char* chain4 = "<GRAPH name=\"chain4\" vertices=\"4\"><EDGE source=\"1\" target=\"2\"/>"
                     "<EDGE source=\"2\" target=\"3\"/><EDGE source=\"3\" target=\"4\"/></GRAPH>";

TEST(SectorCheck, NonConservingBondTermThrows) {
  Files f(chain4, "<HAMILTONIAN name=\"bond pairing\"><BASIS ref=\"fermion\"/>"
          "<BONDTERM source=\"i\" target=\"j\">-fermion_hop(i,j)"
          " + 0.7*(cdag_up(i)*cdag_down(j) + c_down(j)*c_up(i))</BONDTERM></HAMILTONIAN>");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "chain4"; p["MODEL"] = "bond pairing";
  EXPECT_EQ(spectrum(p).size(), 256u);
  std::string e = build_error(sector(p, 1, 1));
  EXPECT_NE(e.find("<BONDTERM>"), std::string::npos) << e;
  EXPECT_NE(e.find("Nup_total"), std::string::npos) << e;
}

TEST(SectorCheck, NonConservingSiteTermThrows) {
  Files f(chain4, "<HAMILTONIAN name=\"site pairing\"><BASIS ref=\"fermion\"/>"
          "<SITETERM site=\"i\">0.3*(cdag_up(i)*cdag_down(i) + c_down(i)*c_up(i))</SITETERM>"
          "<BONDTERM source=\"i\" target=\"j\">-fermion_hop(i,j)</BONDTERM></HAMILTONIAN>");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "chain4"; p["MODEL"] = "site pairing";
  EXPECT_EQ(spectrum(p).size(), 256u);
  EXPECT_NE(build_error(sector(p, 2, 1)).find("<SITETERM>"), std::string::npos);
}

TEST(SectorCheck, NonConservingPlaquetteTermThrowsInMomentumSectors) {
  Files f(ladder_cell, hubbard_model("plaquette pairing", "<PLAQUETTETERM sites=\"i j k l\">"
          "0.25*(cdag_up(i)*cdag_down(k)*cdag_up(l)*c_up(j) + cdag_up(j)*c_up(l)*c_down(k)*c_up(i))"
          "</PLAQUETTETERM>"));
  alps::Parameters p = sector(f.parms(), 3, 3);
  p["LATTICE"] = "periodic ladder"; p["L"] = 3; p["MODEL"] = "plaquette pairing";
  for (const char* k : {"0", "2.0943951023931953"}) {
    p["TOTAL_MOMENTUM"] = k;
    EXPECT_NE(build_error(p).find("<PLAQUETTETERM>"), std::string::npos) << "k=" << k;
  }
}

TEST(SectorCheck, ConservingTermsInMomentumSectorsDoNotThrow) {
  // momentum sectors drop states whose orbit has zero norm; that is not an error
  Files f(ladder_cell, hubbard_model("pair hopping",
          "<PLAQUETTETERM sites=\"i j k l\">" + std::string(pair_hop) + "</PLAQUETTETERM>"));
  alps::Parameters p = sector(f.parms(), 3, 2);
  p["LATTICE"] = "periodic ladder"; p["L"] = 3; p["MODEL"] = "pair hopping";
  p["U"] = 8.0; p["K"] = 0.9;
  for (const char* k : {"0", "2.0943951023931953", "-2.0943951023931953"}) {
    p["TOTAL_MOMENTUM"] = k;
    EXPECT_EQ(build_error(p), "") << "k=" << k;
  }
}

TEST(SectorCheck, MeasuringANonConservingOperatorGivesZero) {
  // operator_matrix (measurements) keeps dropping elements that leave the sector
  Files f(chain4, "<HAMILTONIAN name=\"hop\"><BASIS ref=\"fermion\"/>"
          "<BONDTERM source=\"i\" target=\"j\">-fermion_hop(i,j)</BONDTERM></HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 1, 1);
  p["GRAPH"] = "chain4"; p["MODEL"] = "hop";
  typedef alps::numeric::matrix<double> matrix;
  alps::hamiltonian_matrix<matrix> h(p);
  matrix m = h.operator_matrix<matrix>(alps::SiteOperator("cdag_up(i)*cdag_down(i)", "i"), 1);
  for (std::size_t i = 0; i < h.dimension(); ++i)
    for (std::size_t j = 0; j < h.dimension(); ++j)
      EXPECT_EQ(m(i, j), 0.);
}

TEST(NSiteHamiltonian, TableIndexOverflowThrows) {
  // The local table keys each configuration of the plaquette's sites by one
  // std::size_t; with 32 fermion sites (4^32 = 2^64) it would wrap.
  Files f("", "");
  alps::ModelLibrary lib(f.parms());
  for (int n : {31, 32}) {
    std::vector<std::string> roles;
    for (int r = 0; r < n; ++r)
      roles.push_back("r" + std::to_string(r));
    const std::string last = roles.back();
    alps::NSiteOperator op("cdag_up(r0)*c_up(" + last + ") + cdag_up(" + last + ")*c_up(r0)", roles);
    std::vector<alps::SiteBasisDescriptor<short> > bases(n, lib.get_site_basis("fermion"));
    if (n == 31)
      EXPECT_NO_THROW(op.matrix<double>(bases, f.parms()));
    else
      EXPECT_THROW(op.matrix<double>(bases, f.parms()), std::runtime_error);
  }
}

TEST(NSiteHamiltonian, DenseMatrixSizeOverflowThrows) {
  // 2^33 x 2^33 = 2^66 wraps to 0 in a std::size_t
  const std::size_t d = std::size_t(1) << 33;
  EXPECT_THROW(alps::detail::initialized_matrix(alps::numeric::matrix<double>(), d, d), std::runtime_error);
  EXPECT_EQ(num_rows(alps::detail::initialized_matrix(alps::numeric::matrix<double>(), 3, 3)), 3u);
}

TEST(NSiteHamiltonian, TypedAndDefaultTerms) {
  // Two plaquette types; '#' in the untyped term resolves per type.
  Files f("<GRAPH name=\"g\" vertices=\"3\">"
          "<PLAQUETTE type=\"0\"><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE>"
          "<PLAQUETTE type=\"1\"><SITE id=\"2\"/><SITE id=\"3\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"typed\"><BASIS ref=\"fermion\"/>"
          "<PLAQUETTETERM sites=\"a b\"><PARAMETER name=\"V#\" default=\"0\"/>V#*n(a)*n(b)</PLAQUETTETERM>"
          "</HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 3, 0);
  p["GRAPH"] = "g"; p["MODEL"] = "typed"; p["V0"] = 0.25; p["V1"] = 2.0;
  std::vector<double> e = spectrum(p);
  ASSERT_EQ(e.size(), 1u);
  EXPECT_NEAR(e[0], 2.25, 1e-14);
}

TEST(NSiteHamiltonian, TypedTermKeepsItsParameterDefaults) {
  Files f("<GRAPH name=\"g\" vertices=\"2\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"typed\"><BASIS ref=\"fermion\"/>"
          "<PLAQUETTETERM type=\"0\" sites=\"a b\"><PARAMETER name=\"K\" default=\"0.7\"/>K*n(a)*n(b)</PLAQUETTETERM>"
          "</HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 2, 0);
  p["GRAPH"] = "g"; p["MODEL"] = "typed";
  EXPECT_NEAR(spectrum(p)[0], 0.7, 1e-14);
  p["K"] = 2.0;                           // a given parameter wins over the default
  EXPECT_NEAR(spectrum(p)[0], 2.0, 1e-14);
}

// The splitter replaces each site operator by a sign, so n(a)^2 used to act
// as n(a), and sqrt(n(a)) as n(a).
TEST(NSiteHamiltonian, PowersAndFunctionsOfSiteOperatorsThrow) {
  const char* g = "<GRAPH name=\"g\" vertices=\"2\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE></GRAPH>";
  Files f(g, "<HAMILTONIAN name=\"product\"><BASIS ref=\"fermion\"/>"
             "<PLAQUETTETERM sites=\"a b\">n(a)*n(a)*n(b)</PLAQUETTETERM></HAMILTONIAN>"
             "<HAMILTONIAN name=\"power\"><BASIS ref=\"fermion\"/>"
             "<PLAQUETTETERM sites=\"a b\">n(a)^2*n(b)</PLAQUETTETERM></HAMILTONIAN>"
             "<HAMILTONIAN name=\"function\"><BASIS ref=\"fermion\"/>"
             "<PLAQUETTETERM sites=\"a b\">sqrt(n(a))*n(b)</PLAQUETTETERM></HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 2, 1);   // occupations (2,1) and (1,2)
  p["GRAPH"] = "g"; p["MODEL"] = "product";
  std::vector<double> e = spectrum(p);
  ASSERT_EQ(e.size(), 2u);
  EXPECT_NEAR(e[0], 2., 1e-14);
  EXPECT_NEAR(e[1], 4., 1e-14);
  for (const char* m : {"power", "function"}) {
    p["MODEL"] = m;
    EXPECT_THROW(spectrum(p), std::runtime_error) << m;
  }
}

TEST(NSiteHamiltonian, ApplyOperatorBuildsTheBasis) {
  Files f("<GRAPH name=\"g\" vertices=\"2\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"m\"><BASIS ref=\"fermion\"/>"
          "<PLAQUETTETERM sites=\"a b\">n(a)*n(b)</PLAQUETTETERM></HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 2, 0);   // one state, n = 1 on both sites
  p["GRAPH"] = "g"; p["MODEL"] = "m";
  alps::hamiltonian_matrix<alps::numeric::matrix<double> > h(p);
  alps::plaquette plaq;
  plaq.sites = {0, 1};
  std::vector<double> x(1, 1.), y(1, 0.);
  h.apply_operator(alps::NSiteOperator("n(a)*n(b)", {"a", "b"}), plaq, x, y);
  ASSERT_EQ(h.dimension(), 1u);
  EXPECT_NEAR(y[0], 1., 1e-14);
}

TEST(NSiteHamiltonian, ParametersKeepFullPrecision) {
  // Substituted parameters used to be printed with 6 significant digits.
  Files f("<GRAPH name=\"g\" vertices=\"2\"><EDGE source=\"1\" target=\"2\"/>"
          "<PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"m\"><PARAMETER name=\"t\" default=\"0\"/><PARAMETER name=\"V\" default=\"0\"/>"
          "<BASIS ref=\"fermion\"/><BONDTERM source=\"i\" target=\"j\">-t*fermion_hop(i,j)</BONDTERM>"
          "<PLAQUETTETERM sites=\"a b\">V*n(a)*n(b) + 0.123456789123*n_up(a)</PLAQUETTETERM></HAMILTONIAN>");
  alps::Parameters p = sector(f.parms(), 1, 0);
  p["GRAPH"] = "g"; p["MODEL"] = "m"; p["t"] = 0.987654321987;
  std::vector<double> e = spectrum(p);
  const double d = 0.123456789123, t = 0.987654321987;
  EXPECT_NEAR(e[0], d / 2 - std::sqrt(d * d / 4 + t * t), 1e-13);
}

TEST(NSiteHamiltonian, ApplicationsWithoutSupportRefuse) {
  Files f(ladder_cell, hubbard_model("pair hopping",
          "<PLAQUETTETERM sites=\"i j k l\">" + std::string(pair_hop) + "</PLAQUETTETERM>"));
  alps::Parameters p = f.parms();
  p["LATTICE"] = "open ladder"; p["L"] = 2; p["MODEL"] = "pair hopping";
  alps::graph_helper<> g(p);
  EXPECT_THROW(alps::model_helper<> m(g, p), std::runtime_error);
  EXPECT_NO_THROW(alps::model_helper<> m(g, p, false, true));
}

// Bosons: no Jordan-Wigner sign, local range N = 0..Nmax. A term that would
// take a site above Nmax has that matrix element dropped. The reference below
// builds H densely on all (Nmax+1)^4 occupations of a 4-site ring.
namespace {

const char* boson_ring =
  "<GRAPH name=\"ring\" vertices=\"4\">"
  "<EDGE source=\"1\" target=\"2\"/><EDGE source=\"2\" target=\"3\"/>"
  "<EDGE source=\"3\" target=\"4\"/><EDGE source=\"4\" target=\"1\"/>"
  "<PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/><SITE id=\"3\"/><SITE id=\"4\"/></PLAQUETTE>"
  "</GRAPH>";

std::string boson_model(const std::string& name, const std::string& plaquette_term,
                        const std::string& sites = "i j k l") {
  return "<SITEBASIS name=\"boson\"><PARAMETER name=\"Nmax\" default=\"infinity\"/>"
         "<QUANTUMNUMBER name=\"N\" min=\"0\" max=\"Nmax\"/>"
         "<OPERATOR name=\"bdag\" matrixelement=\"sqrt(N+1)\"><CHANGE quantumnumber=\"N\" change=\"1\"/></OPERATOR>"
         "<OPERATOR name=\"b\" matrixelement=\"sqrt(N)\"><CHANGE quantumnumber=\"N\" change=\"-1\"/></OPERATOR>"
         "<OPERATOR name=\"n\" matrixelement=\"N\"/></SITEBASIS>"
         "<BASIS name=\"boson\"><SITEBASIS ref=\"boson\"/><CONSTRAINT quantumnumber=\"N\" value=\"N_total\"/></BASIS>"
         "<HAMILTONIAN name=\"" + name + "\">"
         "<PARAMETER name=\"t\" default=\"1\"/><PARAMETER name=\"U\" default=\"0\"/>"
         "<PARAMETER name=\"K\" default=\"0\"/><BASIS ref=\"boson\"/>"
         "<SITETERM site=\"i\">U/2*n(i)*(n(i)-1)</SITETERM>"
         "<BONDTERM source=\"i\" target=\"j\">-t*(bdag(i)*b(j)+bdag(j)*b(i))</BONDTERM>"
         "<PLAQUETTETERM sites=\"" + sites + "\">" + plaquette_term + "</PLAQUETTETERM></HAMILTONIAN>";
}

// One product of boson operators: (site, +1 for bdag / -1 for b), leftmost first.
typedef std::vector<std::pair<int, int> > boson_product;

// Dense H of a 4-site ring on occupations 0..nmax, restricted to total N
// = n_total (all states if n_total < 0). Sorted eigenvalues.
std::vector<double> boson_reference(int nmax, int n_total, double t, double u, double k,
                                    const std::vector<boson_product>& plaquette) {
  const int L = 4, d = nmax + 1;
  std::vector<std::vector<int> > states;
  for (int s = 0; s < d * d * d * d; ++s) {
    std::vector<int> occ(L);
    int n = 0;
    for (int r = 0, x = s; r < L; ++r, x /= d)
      n += occ[r] = x % d;
    if (n_total < 0 || n == n_total)
      states.push_back(occ);
  }
  std::vector<std::vector<double> > h(states.size(), std::vector<double>(states.size()));
  auto add = [&](std::size_t col, double coef, const boson_product& ops) {
    std::vector<int> occ = states[col];
    double amp = coef;
    for (std::size_t o = ops.size(); o-- > 0;) {
      int& m = occ[ops[o].first];
      amp *= ops[o].second > 0 ? std::sqrt(m + 1.0) : std::sqrt(double(m));
      m += ops[o].second;
      if (m < 0)
        return;
    }
    for (int m : occ)
      if (m > nmax)  // outside the truncated basis
        return;
    std::size_t row = std::find(states.begin(), states.end(), occ) - states.begin();
    ASSERT_LT(row, states.size()) << "reference term leaves the sector";
    h[row][col] += amp;
  };
  for (std::size_t c = 0; c < states.size(); ++c) {
    for (int s = 0; s < L; ++s)
      h[c][c] += u / 2 * states[c][s] * (states[c][s] - 1);
    for (int s = 0; s < L; ++s) {
      add(c, -t, {{s, 1}, {(s + 1) % L, -1}});
      add(c, -t, {{(s + 1) % L, 1}, {s, -1}});
    }
    for (const boson_product& ops : plaquette)
      add(c, k, ops);
  }
  return eigenvalues(h);
}

alps::Parameters boson_parms(const Files& f, int nmax) {
  alps::Parameters p = f.parms();
  p["GRAPH"] = "ring"; p["MODEL"] = "m";
  p["Nmax"] = nmax; p["t"] = 1.0; p["U"] = 2.0; p["K"] = 0.7;
  return p;
}

} // namespace

TEST(NSiteBoson, ConservingTermMatchesReference) {
  // Correlated pair hopping (i,j) <-> (k,l); conserves N, so no cutoff warning.
  Files f(boson_ring, boson_model("m", "K*(bdag(i)*bdag(j)*b(k)*b(l) + bdag(l)*bdag(k)*b(j)*b(i))"));
  const std::vector<boson_product> term = {{{0, 1}, {1, 1}, {2, -1}, {3, -1}},
                                           {{3, 1}, {2, 1}, {1, -1}, {0, -1}}};
  for (int n_total = 2; n_total <= 4; ++n_total) {
    alps::Parameters p = boson_parms(f, 2);
    p["CONSERVED_QUANTUMNUMBERS"] = "N";
    p["N_total"] = n_total;
    testing::internal::CaptureStderr();
    std::vector<double> e = spectrum(p);
    EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
    std::vector<double> expected = boson_reference(2, n_total, 1.0, 2.0, 0.7, term);
    ASSERT_EQ(e.size(), expected.size()) << "N_total " << n_total;
    for (std::size_t i = 0; i < e.size(); ++i)
      EXPECT_NEAR(e[i], expected[i], 1e-10) << "N_total " << n_total << " level " << i;
  }
}

TEST(NSiteBoson, NonConservingTermIsTruncatedAtNmaxAndWarns) {
  // Pair creation/annihilation on (i,k); the basis is all occupations 0..Nmax.
  Files f(boson_ring, boson_model("m", "K*(bdag(i)*bdag(k) + b(k)*b(i))"));
  const std::vector<boson_product> term = {{{0, 1}, {2, 1}}, {{2, -1}, {0, -1}}};
  for (int nmax = 1; nmax <= 3; ++nmax) {
    alps::Parameters p = boson_parms(f, nmax);
    testing::internal::CaptureStderr();
    std::vector<double> e = spectrum(p);
    std::string warning = testing::internal::GetCapturedStderr();
    EXPECT_NE(warning.find("Warning: a <PLAQUETTETERM> changes the total N"), std::string::npos) << warning;
    EXPECT_NE(warning.find("N <= Nmax = " + std::to_string(nmax)), std::string::npos) << warning;
    std::vector<double> expected = boson_reference(nmax, -1, 1.0, 2.0, 0.7, term);
    ASSERT_EQ(e.size(), expected.size()) << "Nmax " << nmax;
    for (std::size_t i = 0; i < e.size(); ++i)
      EXPECT_NEAR(e[i], expected[i], 1e-10) << "Nmax " << nmax << " level " << i;
  }
  // Conserving N is impossible with this term.
  alps::Parameters p = boson_parms(f, 2);
  p["CONSERVED_QUANTUMNUMBERS"] = "N";
  p["N_total"] = 2;
  EXPECT_THROW(spectrum(p), std::runtime_error);
}

TEST(NSiteBoson, EachOperatorIsTruncated) {
  // bdag on N = Nmax is zero, so b*bdag there is 0, not Nmax+1.
  Files f("<GRAPH name=\"ring\" vertices=\"1\"><PLAQUETTE><SITE id=\"1\"/></PLAQUETTE></GRAPH>",
          boson_model("m", "K*b(i)*bdag(i)", "i"));
  alps::Parameters p = boson_parms(f, 1);
  p["K"] = 1.0;
  testing::internal::CaptureStderr();
  std::vector<double> e = spectrum(p);
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
  ASSERT_EQ(e.size(), 2u);
  EXPECT_NEAR(e[0], 0.0, 1e-14);  // N = 1 = Nmax
  EXPECT_NEAR(e[1], 1.0, 1e-14);  // N = 0: b*bdag = 1
}

TEST(NSiteBoson, NmaxUnsetIsAnError) {
  Files f(boson_ring, boson_model("m", "K*(bdag(i)*bdag(k) + b(k)*b(i))"));
  alps::Parameters p = f.parms();
  p["GRAPH"] = "ring"; p["MODEL"] = "m";
  EXPECT_THROW(spectrum(p), std::runtime_error);
}

// Totals over all sites are half_integer<int>; they used to be
// half_integer<short> (|total| <= 16383) and wrapped silently. Building these
// bases also needs the enumeration to prune every infeasible value of a site,
// not just one.
std::size_t basis_dimension(const alps::Parameters& p) {
  alps::hamiltonian_matrix<alps::numeric::matrix<double> > h(p);
  return h.dimension();
}

TEST(NSiteBoson, TotalsBeyondShortRange) {
  // 4 sites x Nmax 5000: totals up to 20000
  Files f(boson_ring, boson_model("m", "K*n(i)*n(j)"));
  alps::Parameters p = f.parms();
  p["GRAPH"] = "ring"; p["MODEL"] = "m"; p["Nmax"] = 5000;
  p["CONSERVED_QUANTUMNUMBERS"] = "N";
  p["N_total"] = 2;
  EXPECT_EQ(basis_dimension(p), 10u);       // C(5,3)
  p["N_total"] = 19999;
  EXPECT_EQ(basis_dimension(p), 4u);        // one site at 4999
  p["N_total"] = 20000;
  EXPECT_EQ(basis_dimension(p), 1u);        // all at Nmax
  p["N_total"] = 20001;
  EXPECT_EQ(basis_dimension(p), 0u);
}

TEST(NSiteBoson, LargeSiteRangeWarns) {
  // a single site's range cannot be widened: warn
  Files f(boson_ring, boson_model("m", "K*n(i)*n(j)"));
  alps::Parameters p = f.parms();
  p["GRAPH"] = "ring"; p["MODEL"] = "m"; p["Nmax"] = 20000;
  p["CONSERVED_QUANTUMNUMBERS"] = "N"; p["N_total"] = 1;
  testing::internal::CaptureStderr();
  try { basis_dimension(p); } catch (std::exception&) {}
  std::string err = testing::internal::GetCapturedStderr();
  EXPECT_NE(err.find("maximum of quantum number N"), std::string::npos) << err;
}

TEST(NSiteHamiltonian, SparseSectorOnManySitesIsFast) {
  // 31 sites, one or three fermions: the enumeration used to walk ~4^31 states
  std::ostringstream graph;
  graph << "<GRAPH name=\"chain31\" vertices=\"31\">";
  for (int i = 1; i < 31; ++i)
    graph << "<EDGE source=\"" << i << "\" target=\"" << i + 1 << "\"/>";
  graph << "</GRAPH>";
  Files f(graph.str(), "<HAMILTONIAN name=\"hop\"><BASIS ref=\"fermion\"/>"
          "<BONDTERM source=\"i\" target=\"j\">-fermion_hop(i,j)</BONDTERM></HAMILTONIAN>");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "chain31"; p["MODEL"] = "hop";
  EXPECT_EQ(basis_dimension(sector(p, 1, 0)), 31u);
  EXPECT_EQ(basis_dimension(sector(p, 2, 1)), 465u * 31u);
}

TEST(NSiteBoson, NoCutoffWarningForFermions) {
  // Fermion occupations are bounded by Pauli, not by a cutoff parameter.
  Files f("<GRAPH name=\"g\" vertices=\"4\"><PLAQUETTE><SITE id=\"1\"/><SITE id=\"2\"/><SITE id=\"3\"/><SITE id=\"4\"/></PLAQUETTE></GRAPH>",
          "<HAMILTONIAN name=\"pairing\"><BASIS ref=\"fermion\"/><PLAQUETTETERM sites=\"i j k l\">"
          "cdag_up(i)*cdag_down(j)*cdag_up(k)*cdag_down(l) + c_down(l)*c_up(k)*c_down(j)*c_up(i)"
          "</PLAQUETTETERM></HAMILTONIAN>");
  alps::Parameters p = f.parms();
  p["GRAPH"] = "g"; p["MODEL"] = "pairing";
  testing::internal::CaptureStderr();
  spectrum(p);
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
}
