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

/* $Id$ */

#ifndef ALPS_DMFT_EXTERNALSOLVER_H
#define ALPS_DMFT_EXTERNALSOLVER_H

/// @file externalsolver.h
/// @brief declares the external solver
/// @sa ExternalSolver

#include "solver.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <cstdio>

/// @brief An impurity solver calling an external executable to solve the impurity problem
/// 
/// The ExternalSolver receives the path to an executable in the constructor. This executable
/// will be called to solve the impurity problem. The executable should take two arguments:
/// the name of an input file and the name of an output file. Both files have to conform to the
/// XML schema for impurity solvers that is currently being developed
///
/// @todo the XML formats need to be discussed and documented. responsible: Matthias Troyer
  
class ExternalSolver : public ImpuritySolver {
public:
  /// @param executable the path to the executable
  ExternalSolver(const boost::filesystem::path& executable) 
    : exe_(boost::filesystem::complete(executable)) {}

  multiple_vector_type solve(const multiple_vector_type& G0, const alps::Parameters& parms);
    
private:
  boost::filesystem::path exe_;
};

#endif
