/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include <alps/utility/copyright.hpp>
#include <iostream>

#include "dmrg/version.h"

#include <complex>
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mpo_contractor.h"
#include "dmrg/utils/DmrgParameters.h"

#include <boost/program_options.hpp>


/// build with NU1 symmetry
typedef NU1 grp;


template <class Matrix, class SymmGroup>
void run (std::string const& chkp_in, std::string const& chkp_out, unsigned int maxbond, unsigned int maxiter, double tol)
{
    /// Load state
    MPS<Matrix, grp> mps;
    load(chkp_in, mps);
    
    /// Build identity MPO
    MPO<Matrix, grp> mpo(mps.length());
    for (int p=0; p<mps.length(); ++p) {
        Index<grp> const& phys_i = mps[p].site_dim();
        block_matrix<Matrix, grp> op = identity_matrix<Matrix>(phys_i);
        mpo[p] = MPOTensor<Matrix, grp>(1, 1);
        mpo[p].set(0, 0, op, 1.);
    }
    
    DmrgParameters parms;
    parms.set("max_bond_dimension", maxbond);
    parms.set("mpo_compression", "twosite");
    
    double eps;
    
    /// Compressing the MPS
    mpo_contractor<Matrix, grp, storage::nop> compressor(mps, mpo, parms);
    for (int k = 0; k < maxiter; ++k) {
        std::pair<double,double> sweep_eps = compressor.sweep(k);
        double rel_error = std::abs( (sweep_eps.first - sweep_eps.second) / sweep_eps.second );
        if (rel_error < tol) {
            eps = maquis::real(norm(mps)) + sweep_eps.second;
            break;
        }
        if (k == maxiter-1) {
            std::cerr << "Accuracy of variational compression reached only " << rel_error << " (maxiter="<< maxiter << ", tol=" << tol << ")" << std::endl;
            exit(1); // TODO: nicer to exit with exception and catch it
        }
    }
    compressor.finalize();
    MPS<Matrix, grp> const& mps_compressed = compressor.get_current_mps();
    
    /// Compute resulting max bond dimension
    std::size_t bond_dim = 0;
    for(int p=0; p<mps_compressed.length(); ++p) {
        bond_dim = std::max(mps_compressed[p].col_dim().sum_of_sizes(), bond_dim);
    }
    
    /// Report outcome
    std::cout << "bond_dim = " << bond_dim << std::endl;
    std::cout << "eps      = " << eps << std::endl;
    
    if (!chkp_out.empty()) {
        if (!boost::filesystem::exists(chkp_out))
            boost::filesystem::create_directory(chkp_out);
        save(chkp_out, mps_compressed);
        std::cout << "chkp     = " << chkp_out << std::endl;
    }
}

int main(int argc, char ** argv)
{
    try {
        std::cout << "ALPS/MPS version " DMRG_VERSION_STRING " (2013-2014)\n"
                  << "  Density Matrix Renormalization Group algorithm\n"
                  << "  available from http://alps.comp-phys.org/\n"
                  << "  copyright (c) 2013 Institute for Theoretical Physics, ETH Zurich\n"
                  << "  copyright (c) 2010-2011 by Bela Bauer\n"
                  << "  copyright (c) 2011-2013 by Michele Dolfi\n"
                  << "  for details see the publication: \n"
                  << "  M. Dolfi et al., Computer Physics Communications 185, 3430 (2014).\n"
                  << "                   doi: 10.1016/j.cpc.2014.08.019\n"
                  << std::endl;
        alps::print_copyright(std::cout);
        
        /// parse options
        std::string chkp_in, chkp_out;
        unsigned int maxbond;
        unsigned int maxiter;
        double tol;
        
        namespace po = boost::program_options;
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version", "print program version")
        ("complex", "use complex numbers")
        ("maxbond", po::value<unsigned int>(&maxbond)->required(), "maximum bond dimension")
        ("maxiter", po::value<unsigned int>(&maxiter)->default_value(4), "maximum number of iteration")
        ("tol", po::value<double>(&tol)->default_value(1e-6), "relative tolerance on eps")
        ("input", po::value<std::string>(&chkp_in)->required(), "path to chkp of the input")
        ("output", po::value<std::string>(&chkp_out), "path to chkp for the output");
        po::positional_options_description p;
        p.add("input", 1);
        p.add("output", 1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        if (vm.count("version")) {
            std::cout << alps::version_string() << std::endl;
            std::cout << DMRG_VERSION_STRING << std::endl;
            return 0;
        }
        
        po::notify(vm);
        
        /// compute
        if (vm.count("complex")) run<cmatrix, grp>(chkp_in, chkp_out, maxbond, maxiter, tol);
        else                     run<matrix,  grp>(chkp_in, chkp_out, maxbond, maxiter, tol);
        
    } catch (std::exception & e) {
        std::cerr << "Exception thrown:" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}
