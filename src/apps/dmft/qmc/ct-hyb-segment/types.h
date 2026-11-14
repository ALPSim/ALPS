 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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

//#define uint unsigned int 
#include<sys/types.h>
#include<complex>
/// @file types.h
/// @brief definition for basic types like vectors that are used throughout the solvers.
#ifndef TYPES_H_
#define TYPES_H_
//includes
#include <vector>
#include <complex>
/// the array (vector) type to store frequence dependent Green's functions and Weiss fields
typedef std::vector<double> vector_type;
/// two arrays are needed for a single site problem: one for up and one for down spins
typedef std::pair<vector_type,vector_type> multiple_vector_type;
/// variant if you want complex vectors
typedef std::vector<std::complex<double> > complex_vector_type;
/// variant for two complex vectors.
typedef std::pair<complex_vector_type,complex_vector_type> multiple_complex_vector_type;

///dense matrix
//typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
///complex matrix
//typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_matrix;

///same as complex_vector_type
typedef std::vector<std::complex<double> > complex_vector;
///same as vector_type
typedef std::vector<double> double_vector;
///vector of ints
typedef std::vector<int> int_vector;

///enum for spin up and spin down
enum  {up=0, down=1} ;
///addressing type for site indices (cluster)
typedef unsigned int site_t;
///addressing type for spin indices
typedef unsigned int spin_t;
///type of imaginary time values
typedef double itime_t;
///addressing type of imaginary time indices (discretized)
typedef unsigned int itime_index_t;
///addressing type of matsubara frequency
typedef unsigned int frequency_t;

///typedef to easily distinguish imaginary time and frequency (multiple) vectors
typedef  multiple_vector_type itime_multiple_vector_t;
///typedef to easily distinguish imaginary time and frequency (multiple) vectors
typedef  multiple_complex_vector_type freq_multiple_vector_t;


std::ostream &operator<<(std::ostream &os, const multiple_vector_type &v);
std::ostream &operator<<(std::ostream &os, const multiple_complex_vector_type &v);
std::ostream &operator<<(std::ostream &os, const double_vector &v);
std::ostream &operator<<(std::ostream &os, const complex_vector &v);
void print_imag_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta);
void print_real_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta);
  
void print_green_itime(std::ostream &os, const multiple_vector_type &v, const double beta);

struct hifreq_moments{
  std::vector<std::vector<std::vector<double> > >c0;
  std::vector<std::vector<std::vector<double> > >c1;
  std::vector<std::vector<std::vector<double> > >c2;
  std::vector<std::vector<std::vector<double> > >c3;
};
#endif /*TYPES_H_*/
