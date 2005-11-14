/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
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

//definition for basic types like vectors that are used throughout the solvers.
#ifndef TYPES_H_
#define TYPES_H_

//includes
#include <vector>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

/// the array (vector) type to store frequence dependent Green's functions and Weiss fields
typedef std::vector<double> vector_type;
/// two arrays are needed for a single site problem: one for up and one for down spins
typedef std::pair<vector_type,vector_type> multiple_vector_type;
/// variant if you want complex vectors
typedef std::vector<std::complex<double> > complex_vector_type;
/// variant for two complex vectors.
typedef std::pair<complex_vector_type,complex_vector_type> multiple_complex_vector_type;

///dense matrix
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;



#endif /*TYPES_H_*/
