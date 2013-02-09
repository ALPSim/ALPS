/***************************************************************************
 * $Id: jacobi.cpp,v 1.9 2004/06/29 08:31:02 troyer Exp $
 *
 * An example of the Lanczos method for the calculation of n lowest eigenvalues.
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
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
 *
 **************************************************************************/

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/std.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>
#include <boost/numeric/bindings/upper.hpp> 

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/jacobi.h>
#include <ietl/iteration.h>

#include <boost/random.hpp>
#include <boost/limits.hpp>

#include <cmath>
#include <limits>

typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> Matrix; 
typedef boost::numeric::ublas::vector<double> Vector;

int main() {
    // Creation of an example matrix:
    int N = 500;
    Matrix mat(N, N);
    int n = 1;
    for(int i=0;i<N;i++)
        for(int j=0;j<=i;j++) {
            mat(i,j) = n++;
            // mat(i,j) = drand48(); // n++;
            // if (i == j)
            //     mat(i,i) += 10;
            mat(j,i) = mat(i,j);
        }
    
    typedef ietl::vectorspace<Vector> Vecspace;
    typedef boost::lagged_fibonacci607 Gen;  
    
    Vecspace vec(N);
    Gen mygen;
    
    ietl::jcd_simple_solver<Matrix, Vecspace> jcd_ss(mat, vec);
    ietl::jcd_gmres_solver<Matrix, Vecspace> jcd_gmres(mat, vec);
    ietl::jacobi_davidson<Matrix, Vecspace> jd(mat, vec, ietl::Smallest);
    
    ietl::basic_iteration<double> iter(200, 1e-8, 1e-8);
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(mygen, jcd_gmres, iter);
    std::cout << "Jacobi-Davidson says: " << r0.first << std::endl;
    std::cout << "It used " << iter.iterations() << " iterations." << std::endl;

    Vector v2 = new_vector(vec);
    ietl::mult(mat, r0.second, v2);
    std::cout << "Output vector norm: " << ietl::two_norm(r0.second) << std::endl;
    std::cout << "Checked eigenvalue: " << ietl::dot(r0.second, v2) << std::endl;
    
    {
        std::cout << "LAPACK says: ";
        Vector evals(N);
        boost::numeric::bindings::lapack::heev('N', boost::numeric::bindings::upper(mat), evals);
    
        std::copy(evals.begin(), evals.begin()+3, std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;
    }
    
    return 0;
}
