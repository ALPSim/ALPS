/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef COMPRESSION_H
#define COMPRESSION_H

#include "dmrg/mp_tensors/mps.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"

struct compression {
    template<class Matrix, class SymmGroup>
    static truncation_results
    replace_two_sites_l2r(MPS<Matrix, SymmGroup> & mps,
                          std::size_t Mmax, double cutoff,
                          block_matrix<Matrix, SymmGroup> const & t,
                          std::size_t p)
    {
        block_matrix<Matrix, SymmGroup> u, v;
        
        typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
        
        block_matrix<dmt, SymmGroup> s;
        
        truncation_results trunc = svd_truncate(t, u, v, s,
                                                cutoff, Mmax, true);
        
        mps[p].replace_left_paired(u, Lnorm);
        
        gemm(s, v, u);
        mps[p+1].replace_right_paired(u);
        
        return trunc;
    }

    template<class Matrix, class SymmGroup>
    static truncation_results
    replace_two_sites_r2l(MPS<Matrix, SymmGroup> & mps,
                          std::size_t Mmax, double cutoff,
                          block_matrix<Matrix, SymmGroup> const & t,
                          std::size_t p)
    {
        block_matrix<Matrix, SymmGroup> u, v;
        
        typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
        
        block_matrix<dmt, SymmGroup> s;
        
        truncation_results trunc = svd_truncate(t, u, v, s,
                                                cutoff, Mmax, true);
        
        mps[p+1].replace_right_paired(v, Rnorm);
        
        gemm(u, s, v);
        mps[p].replace_left_paired(v);
        
        return trunc;
    }
    
    template<class Matrix, class SymmGroup>
    static void compress_two_sites(MPS<Matrix, SymmGroup> & mps,
                                   std::size_t Mmax, double cutoff,
                                   std::size_t p)
    {
        block_matrix<Matrix, SymmGroup> t;
        
        mps[p].make_left_paired();
        mps[p+1].make_right_paired();
        
        gemm(mps[p].data(), mps[p+1].data(), t);
        
        replace_two_sites_l2r(mps, Mmax, cutoff, t, p);
    }
        
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static l2r_compress(MPS<Matrix, SymmGroup> mps,
                        std::size_t Mmax, double cutoff,
                        bool verbose = false)
    {
        std::size_t L = mps.length();
        std::vector<double> ret;
        
        block_matrix<Matrix, SymmGroup> t;
        
        mps.canonize(1);
        
        if (verbose) maquis::cout << "Compressing @ ";
        for (std::size_t p = 1; p < L; ++p)
        {
            if (verbose) {
                maquis::cout << p << " ";
                maquis::cout.flush();
            }
            
            compress_two_sites(mps, Mmax, cutoff, p-1);
            
            t = mps[p].normalize_left(DefaultSolver());
            
            if (p+1 < L)
                mps[p+1].multiply_from_left(t);
            else
                maquis::cout << "Norm reduction: " << trace(t) << std::endl;
        }
        
        return mps;
    }
    
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static r2l_compress(MPS<Matrix, SymmGroup> mps,
                        std::size_t Mmax, double cutoff,
                        bool verbose = false)
    {
        std::size_t L = mps.length();
        std::vector<double> ret;
        
        block_matrix<Matrix, SymmGroup> t;
        
        mps.canonize(L-1);
        
        if (verbose) maquis::cout << "Compressing @ ";
        for (std::size_t p = L-1; p > 0; --p)
        {
            if (verbose) {
                maquis::cout << p << " ";
                maquis::cout.flush();
            }
            
            compress_two_sites(mps, Mmax, cutoff, p-1);
            
            t = mps[p-1].normalize_right(DefaultSolver());
            
            if (p > 1)
                mps[p-2].multiply_from_right(t);
            else
                maquis::cout << "Norm reduction: " << trace(t) << std::endl;
        }
        
        return mps;
    }
};

#endif
