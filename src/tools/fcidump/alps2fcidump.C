/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// alps2fcidump: write an ALPS fermion model as an FCIDUMP file.
//
// The work is done by alps::fcidump::export_fcidump in
// <alps/model/fcidump.hpp>, which also serves the FCIDUMP parameter of
// fulldiag and sparsediag. See that header for how the integrals are read
// off ALPS's matrix and checked.

#include <alps/model/fcidump.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void usage(std::ostream& os)
{
  os << "Usage: alps2fcidump [options] PARAMETER_FILE [FCIDUMP]\n\n"
        "Write the Hamiltonian of an ALPS fermion model (BASIS \"fermion\": one\n"
        "orbital per site; site, bond and n-site terms) as an FCIDUMP file.\n"
        "PARAMETER_FILE is an ALPS parameter file (LATTICE, MODEL, couplings, and\n"
        "LATTICE_LIBRARY / MODEL_LIBRARY if needed); the first parameter set is\n"
        "used. NELEC and MS2 come from Nup_total and Ndown_total, or N_total\n"
        "and Sz_total, if set; otherwise half filling is written.\n"
        "FCIDUMP defaults to 'FCIDUMP'.\n\n"
        "Options:\n"
        "  --check N    compare ALPS's matrix with the integrals in every sector\n"
        "               with at most N electrons (default 4, minimum 2)\n"
        "  --tol X      largest allowed difference (default 1e-10)\n"
        "  --cutoff X   omit integrals with |value| <= X (default 1e-14)\n"
        "  --help       show this message\n";
}

}  // namespace

int main(int argc, char** argv)
{
  try {
    int check = 4;
    double tol = 1e-10, cutoff = 1e-14;
    std::vector<std::string> files;
    for (int k = 1; k < argc; ++k) {
      std::string a = argv[k];
      if (a == "--help" || a == "-h") { usage(std::cout); return 0; }
      else if (a == "--check" && k + 1 < argc) check = boost::lexical_cast<int>(argv[++k]);
      else if (a == "--tol" && k + 1 < argc) tol = boost::lexical_cast<double>(argv[++k]);
      else if (a == "--cutoff" && k + 1 < argc) cutoff = boost::lexical_cast<double>(argv[++k]);
      else if (!a.empty() && a[0] == '-') { usage(std::cerr); return 1; }
      else files.push_back(a);
    }
    if (files.empty() || files.size() > 2) { usage(std::cerr); return 1; }
    std::string target = files.size() > 1 ? files[1] : "FCIDUMP";

    alps::ParameterList list;
    {
      std::ifstream in(files[0].c_str());
      if (!in) throw std::runtime_error("cannot read " + files[0]);
      in >> list;
    }
    alps::Parameters p;
    if (list.empty()) {
      std::ifstream in(files[0].c_str());
      in >> p;
    } else {
      p = list[0];
      if (list.size() > 1)
        std::cerr << "Note: " << list.size() << " parameter sets; using the first.\n";
    }
    if (p.defined("Nup_total") != p.defined("Ndown_total"))
      throw std::runtime_error(p.defined("Nup_total")
          ? "Nup_total is set without Ndown_total; set both, or neither for half filling"
          : "Ndown_total is set without Nup_total; set both, or neither for half filling");
    alps::fcidump::summary s = alps::fcidump::export_fcidump(p, target, check, tol, cutoff);
    std::cout << "Wrote " << target << ": " << s.orbitals << " orbitals, NELEC=" << s.nelec
              << ", MS2=" << s.ms2 << ", E_core=" << std::setprecision(12) << s.ecore
              << "\nChecked all sectors with up to " << s.checked
              << " electrons against ALPS: max |difference| = " << std::setprecision(3)
              << s.worst << "\n";
    return 0;
  } catch (std::exception& e) {
    std::cerr << "alps2fcidump: " << e.what() << "\n";
    return 1;
  }
}
