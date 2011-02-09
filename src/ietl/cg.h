/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
 *****************************************************************************/
 
#ifndef IETL_CG_H
#define IETL_CG_H
 
#include <vector>
#include <iostream>
 
// Note: does not raise any exceptions on failure!
 
namespace ietl
{
    template<class Vector, class Matrix>
    Vector ietl_cg(Matrix const & A, Vector const & b, Vector const & x0, std::size_t max_iter = 100, double abs_tol = 1e-6, bool verbose = false)
    {   
        std::vector<double> rho(max_iter);
    
        // I assume cheap default-constructibility for now
        std::vector<Vector> p(max_iter), q(max_iter), x(max_iter), r(max_iter);
    
        mult(A, x0, x[0]);
        r[0] = b - x[0];
        
        x[0] = x0;
    
        for (std::size_t k = 1; k < max_iter; ++k)
        {
            rho[k-1] = two_norm(r[k-1]);
            rho[k-1] *= rho[k-1];
            if (rho[k-1] < abs_tol)
                return x[k-1];
            
            if (k == 1)
                p[1] = r[0];
            else {
                double beta = rho[k-1]/rho[k-2];
                p[k] = r[k-1] + beta*p[k-1];
            }
        
            mult(A, p[k], q[k]);
            double alpha = rho[k-1]/dot(p[k], q[k]);
            x[k] = x[k-1] + alpha*p[k];
            r[k] = r[k-1] - alpha*q[k];
        
            Vector diff = x[k] - x[k-1];
            double resid = two_norm(diff);
            if (verbose)
                std::cout << "Iteration " << k << ", resid = " << resid << std::endl;
            if (resid < abs_tol)
                return x[k];
        }
        
        return x[max_iter-1];
    }
}
 
#endif

