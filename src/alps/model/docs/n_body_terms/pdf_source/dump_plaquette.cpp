// Scratch: print ALPS's basis, plaquette table and H for the 2x2 plaquette walkthrough.
#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>
#include <iostream>
#include <iomanip>

typedef alps::numeric::matrix<double> M;

int main(int argc, char** argv) {
  alps::Parameters p;
  p["LATTICE_LIBRARY"] = argv[1]; p["MODEL_LIBRARY"] = argv[2];
  p["GRAPH"] = "2x2 plaquette"; p["MODEL"] = "Hubbard with pair hopping";
  p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown"; p["Nup_total"] = 1; p["Ndown_total"] = 1;
  p["t"] = 1.0; p["U"] = 4.0; p["K"] = argv[3];
  alps::hamiltonian_matrix<M> h(p);
  const M& m = h.matrix();
  if (argc > 4) {
    alps::graph_helper<> g(p);
    alps::model_helper<> mh(g, p, false, true);
    alps::site_basis<short> sb(mh.model().basis().site_basis(0));
    std::cout << "LOCAL";
    for (std::size_t k = 0; k < sb.size(); ++k) {
      std::cout << " " << k << ":";
      for (std::size_t q = 0; q < sb[k].size(); ++q) std::cout << sb[k][q] << (q + 1 < sb[k].size() ? "," : "");
    }
    std::cout << "\n";
    const alps::plaquette& pl = g.plaquettes()[0];
    std::cout << "PLAQ";
    for (auto s : pl.sites) std::cout << " " << s;
    std::cout << "\n";
    alps::NSiteOperator term = mh.model().plaquette_term(pl.type);
    std::cout << "TERM " << term.term() << "\n";
    alps::nsite_matrix<double> t = h.local_matrix(term, pl);
    std::size_t nnz = 0;
    std::vector<std::size_t> in(4), out(4);
    for (std::size_t i = 0; i < t.size(); ++i)
      for (auto& e : t.rows[i]) {
        ++nnz; t.decode(i, in); t.decode(e.first, out);
        std::cout << "ENTRY " << in[0] << in[1] << in[2] << in[3] << " " << out[0] << out[1] << out[2] << out[3]
                  << " " << std::setprecision(17) << e.second << "\n";
      }
    std::cout << "TABLE rows " << t.size() << " nnz " << nnz << "\n";
    const auto& st = h.states_vector();
    for (std::size_t i = 0; i < st.size(); ++i) {
      std::cout << "STATE " << i;
      for (std::size_t s = 0; s < st[i].size(); ++s) std::cout << " " << st[i][s];
      std::cout << "\n";
    }
  }
  std::cout << "DIM " << h.dimension() << "\n";
  for (std::size_t i = 0; i < h.dimension(); ++i) {
    std::cout << "ROW";
    for (std::size_t j = 0; j < h.dimension(); ++j) std::cout << " " << std::setprecision(17) << m(i, j);
    std::cout << "\n";
  }
}
