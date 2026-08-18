/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
 
#ifndef IETL_CG_H
#define IETL_CG_H
 
#include <vector>
#include <iostream>
 
// Note: does not raise any exceptions on failure!
 
namespace ietl
{
    class ietl_cg
    {
    private:
        std::size_t max_iter;
        bool verbose;
        
    public:
        ietl_cg(std::size_t max_iter = 100, bool verbose = false)
        : max_iter(max_iter)
        , verbose(verbose) { }
        
        template<class Vector, class Matrix>
        Vector operator()(Matrix const & A, Vector const & b, Vector const & x0, double abs_tol = 1e-6)
        {   
            std::vector<double> rho(2);
        
            // I assume cheap default-constructibility for now
            Vector x, p, q;
        
            mult(A, x0, x);
            Vector r = b - x;
            
            x = x0;
        
            for (std::size_t k = 1; k < max_iter; ++k)
            {
                rho[0] = rho[1];
                rho[1] = two_norm(r);
                rho[1] *= rho[1];
                if (rho[1] < abs_tol)
                    return x;
                
                if (k == 1)
                    p = r;
                else {
                    double beta = rho[1]/rho[0];
                    p = r + beta*p;
                }
            
                mult(A, p, q);
                double alpha = rho[1]/dot(p, q);
                x += alpha*p;
                r -= alpha*q;
            
                Vector diff = alpha*p;
                double resid = two_norm(diff);
                if (verbose)
                    std::cout << "Iteration " << k << ", resid = " << resid << std::endl;
                if (resid < abs_tol)
                    return x;

            }
            
            return x;
        }
    };
}
 
#endif

