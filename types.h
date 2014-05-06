/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*                       Emanuel Gull <gullc@comp-phys.org>
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

#include<complex>
#include <vector>
#include <iostream>

/// @file types.h
/// @brief definition for basic types like vectors that are used throughout the solvers.
#ifndef TYPES_H_
#define TYPES_H_
//includes

///enum for spin up and spin down
enum  {up=0, down=1} ;
///enum for green's function shape
enum shape_t {diagonal, blockdiagonal, nondiagonal};

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

struct hifreq_moments{
  std::vector<std::vector<std::vector<double> > >c0;
  std::vector<std::vector<std::vector<double> > >c1;
  std::vector<std::vector<std::vector<double> > >c2;
  std::vector<std::vector<std::vector<double> > >c3;
};

namespace alps{class Parameters;}

#endif /*TYPES_H_*/
