/***************************************************************************
 * $Id: jacobi.cpp,v 1.9 2004/06/29 08:31:02 troyer Exp $
 *
 * An example of the Lanczos method for the calculation of n lowest eigenvalues.
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
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
 *
 **************************************************************************/

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <boost/timer.hpp>

#include <cassert>

namespace ublas = boost::numeric::ublas;

typedef ublas::symmetric_matrix<double, boost::numeric::ublas::lower> matrix_t;
typedef ublas::banded_matrix<double> b_matrix_t; 
typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double> matrix_type;
typedef boost::lagged_fibonacci44497 gen_t;
typedef boost::normal_distribution<double> dist_t;
// jacobi preconditioner as an example
// the preconditioner must provide a function
// ietl::mult(A, x, b) with b ~= A x
//[ preconditioner
    template <class MATRIX, class SCALAR>
    class jacobi_prec {
        public:
            jacobi_prec( const MATRIX& a, SCALAR l )
             : A_(a), lambda_(l) {}

            template <class VECTOR>
            void mult( const VECTOR & x, VECTOR & y ) const
            {
                assert(A_.size2() == x.size());
                if(y.size() != x.size())
                    y = VECTOR(x.size());

                for(unsigned i = 0; i < x.size(); ++i)
                    y(i) = x(i) / ( A_(i,i) - lambda_);
            }
        private:
            const MATRIX& A_;
            SCALAR lambda_;
    };

namespace ietl {
    template <class MATRIX, class SCALAR, class VECTOR>
    void mult( jacobi_prec<MATRIX,SCALAR> const & K, const VECTOR & x, VECTOR & y )
    {
        K.mult(x, y);
    }
}
//] preconditioner
#include <ietl/jd.h>

int main () {

    int N = 1000;

    // generating a diagonal dominant random matrix A
    dist_t dist(0./*mean*/,1./*var*/);
    gen_t engine;
    boost::variate_generator<gen_t&,dist_t> rng(engine,dist);

    matrix_t A(N, N);

    for(int i=0; i<N ; i++)
        for(int j=0; j < i; j++)
            A(i,j) = rng();

    for(int i=0; i<N; i++)
        A(i,i) = rng() * N;

    //std::cout << "Matrix: " << A << std::endl;

    ietl::vectorspace<vector_t> vs(N);
    gen_t gen;  
    int max_iter = 10*N;
    double rel_tol = 1e-8;
    double abs_tol = 1e-8;
    int m_min = 10, m_max = 20;

    ietl::jd_iteration<double> iter(max_iter, m_min, m_max, rel_tol, abs_tol);
    ietl::jd_iteration<double> iter3(max_iter, m_min, m_max, rel_tol, abs_tol);
    ietl::jd_iteration<double> iter2(max_iter, m_min, m_max, 0.1, abs_tol);

    std::cout.precision(10);  

    ietl::jd<matrix_t, ietl::vectorspace<vector_t> > jd_test  (A, vs);

    int k = 3;

    std::cout << "solve without preconditioning...";
    std::cout.precision(10);
    std::cout.flush();
    boost::timer clock;

    //correction equation solver
    ietl::ietl_gmres solver;

    //search k lowest eigenvalue
    try{
        jd_test.eigensystem(iter, gen, k, solver);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    }
    std::cout << " done. \n\t time: "<< clock.elapsed() << " \t iterations: " << iter.iterations() << "\n";

    jd_test.reset();

    std::cout << "find approximate eigenvalue to create jacobi preconditioner... ";
    std::cout.flush();
    try{
        jd_test.eigensystem(iter2, gen, 1, solver);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    }

    double lambda = jd_test.eigenvalue(0);

    std::cout <<"done. \napprox eigenvalue: "<<lambda<<"\ncreating preconditioner... done.\n";

    //[ funccall
    // create jacobi preconditioner
    jacobi_prec<matrix_t,double> K(A, lambda);

    jd_test.reset();
    // correction equation solver steps are more expensive with preconditioning,
    // but fewer are needed to have 'good' convergence
    ietl::ietl_gmres solver2(3);

    std::cout << "solve with jacobi preconditioning...";
    std::cout.flush();
    clock.restart();
    try{
        jd_test.eigensystem(iter3, gen, k, K, solver2);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    }
    //] funccall

    std::cout << "done. \n time: "<< clock.elapsed() << " \t iterations: " << iter3.iterations() << "\n";

    for(int i = 0; i < jd_test.eigenvalues().size(); ++i)
        std::cout <<"eigenvalue #"<< i <<"\t" << jd_test.eigenvalue(i) <<"\n";/*<<jd_test.eigenvector(i)*///<<"\n\n";

    return 0;
}
