/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
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

/* $Id: externalsolver.h 157 2006-04-04 15:11:22Z gullc $ */

#ifndef ALPS_DMFT_EXTERNALSOLVER_H
#define ALPS_DMFT_EXTERNALSOLVER_H

/// @file externalsolver.h
/// @brief declares the external solver
/// @sa ExternalSolver

#include "solver.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

/// @brief An impurity solver calling an external executable to solve the 
///        impurity problem
/// 
/// The ExternalSolver receives the path to an executable in the constructor. 
/// This executable will be called to solve the impurity problem. The executable
/// should take two arguments: the name of an input file and the name of an 
/// output file. Both files have to conform to the XML schema for impurity 
/// solvers that is currently being developed
///

class ExternalSolver 
 : public ImpuritySolver
 , public MatsubaraImpuritySolver 
{
public:
    /// @param executable the path to the executable
    ExternalSolver(const boost::filesystem::path& executable) ;
  
    ImpuritySolver::result_type solve(
              const itime_green_function_t& G0
            , const alps::Parameters& parms);
    
    MatsubaraImpuritySolver::result_type solve_omega(
              const matsubara_green_function_t& G0_omega
            , const alps::Parameters& parms );
    private:
    /// call the executable
    void call(std::string const& infile, std::string const& outfile);
      
    ///path to the solver executable
    boost::filesystem::path exe_;
};



#endif
