/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Bela Bauer <troyer@comp-phys.org>
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

#include <iostream>
#include <complex>

#include <boost/random.hpp>
#include <boost/limits.hpp>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/geev.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

#include <ietl/interface/ublas.h>
#include <ietl/simple_arnoldi.h>
#include <ietl/vectorspace.h>

inline bool cmp(std::complex<double> a, std::complex<double> b)
{
    return real(a) > real(b);
}

template<class R> double Trand(R& gen, double) { return gen(); }
template<class R> std::complex<double> Trand(R& gen, std::complex<double>) { return std::complex<double>(gen(), gen()); }

template<class T>
void arnoldi1() {
    typedef boost::numeric::ublas::matrix<T, boost::numeric::ublas::column_major> Matrix; 
    typedef boost::numeric::ublas::vector<T> Vector;

    boost::lagged_fibonacci607 gen;

    int N = 10;
    Matrix mat(N, N);
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            mat(i,j) = Trand(gen, T());

    ietl::vectorspace<Vector> vs(N);

    unsigned int nvals = 5;
    ietl::arnoldi_iteration<double> iter(50, nvals, 1e-8, 1e-8);

    ietl::simple_arnoldi<Matrix, ietl::vectorspace<Vector>, boost::lagged_fibonacci607> arni(mat, vs, gen);
    arni.calculate_eigenvalues(iter, true);

    for (unsigned int i = 0; i < nvals; ++i)
       std::cout << "Eigenvalue " << i << ": " << arni.get_eigenvalue(i) << std::endl;

    // Check with dense LAPACK
    {
        Matrix vecs(N,N);
        std::vector<std::complex<double> > evals(N);
        ietl::detail::arnoldi_geev(mat, evals, vecs, T());

        std::sort(evals.begin(), evals.end(), cmp);
        std::copy(evals.begin(), evals.end(), std::ostream_iterator<std::complex<double> >(std::cout, " "));
        std::cout << std::endl;
    }
}
