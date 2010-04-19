/***************************************************************************
 * $Id: simple_arnoldi.h,v 1.34 2004/06/29 08:31:02 troyer Exp $
 *
 * Copyright (C) 2010 by Bela Bauer <bauerb@phys.ethz.ch>
 *
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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
 **************************************************************************/

#ifndef IETL_ARNOLDI_H
#define IETL_ARNOLDI_H

#include <ietl/complex.h>
#include <ietl/iteration.h>
#include <ietl/vectorspace.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>

#include <exception>
#include <stdexcept>
#include <climits>
#include <vector>
#include <iostream>

namespace ietl {
    namespace detail {
        inline bool cmp(std::complex<double> a, std::complex<double> b)
        {
            return ietl::real(a) > ietl::real(b);
        }
    }
    
    template<class T>
    class arnoldi_iteration : public basic_iteration<T>
    {
    public:
        arnoldi_iteration(unsigned int max_iter_, unsigned int desired_eigenvalues__, T reltol_, T abstol_)
            : basic_iteration<T>(max_iter_, reltol_, abstol_)
            , desired_eigenvalues_(desired_eigenvalues__)
        { }
        
        unsigned int desired_eigenvalues() { return desired_eigenvalues_; }
        
    protected:
        unsigned int desired_eigenvalues_;
    };
    
    template<class Matrix, class VS, class Gen>
    class simple_arnoldi
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
        
        simple_arnoldi(Matrix &mat_, VS &vs_, Gen &rng_)
            : mat(mat_)
            , vs(vs_)
            , rng(rng_)
        { }
        
        template<class Iter>
        void calculate_eigenvalues(Iter &iter, bool verbose = false)
        {
            std::vector<vector_type> vectors;
            
            vector_type w = new_vector(vs);
            generate(w, rng);
            project(w, vs);
            w /= two_norm(w);
            
            h_matrix_type H(1,1);
            
            magnitude_type normw, resid;
            
            vectors.push_back(w);
            
            unsigned int j = 0;
            do {
                ietl::mult(mat, vectors[j], w);
                for (int i = 0; i <= j; ++i) {
                    H(i,j) = ietl::dot(w, vectors[i]);
                    w -= H(i,j)*vectors[i];
                }

                normw = ietl::two_norm(w);
                
                // check convergence
                if (j > iter.desired_eigenvalues()) {
                    evals.resize(H.size1());
                    h_matrix_type evecs(H.size1(), H.size2()), H2 = H; // keep a backup because geev destroys the matrix
                    boost::numeric::bindings::lapack::geev(H2, evals, static_cast<h_matrix_type*>(NULL), &evecs,
                        boost::numeric::bindings::lapack::optimal_workspace());
                    double resid = 0;
                    for (int k = 0; k < iter.desired_eigenvalues(); ++k)
                        resid += std::abs(evecs(evecs.size2()-1, k))*normw;
                    if (verbose)
                        std::cout << "Arnoldi iteration " << j << ": residual = " << resid << std::endl;
                    if (iter.finished(resid, abs(evals[iter.desired_eigenvalues()-1]))) {
                        std::sort(evals.begin(), evals.end(), detail::cmp);
                        break;
                    }
                }
                
                H.resize(j+2, j+2, true);
                for (int k = 0; k < j+2; ++k) {
                    H(j+1, k) = 0;
                    H(k, j+1) = 0;
                }
                
                H(j+1,j) = normw;
                w /= normw;
                vectors.push_back(w);
                ++j;
                ++iter;
            } while (true);
        }
        
        std::complex<double> get_eigenvalue(int i)
        {
            return evals[i];
        }
        
    protected:    
        typedef boost::numeric::ublas::matrix<scalar_type, boost::numeric::ublas::column_major> h_matrix_type;
        
        Matrix& mat;
        VS &vs;
        Gen &rng;
        
        std::vector< std::complex<double> > evals;
    };
    
}

#endif
