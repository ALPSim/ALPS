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
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
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
#include <cstdio>

/// @brief An impurity solver calling an external executable to solve the 
///        impurity problem
/// 
/// The ExternalSolver receives the path to an executable in the constructor. 
/// This executable will be called to solve the impurity problem. The executable
/// should take two arguments: the name of an input file and the name of an 
/// output file. Both files have to conform to the XML schema for impurity 
/// solvers that is currently being developed
///
/// @todo the XML formats need to be discussed and documented. 
///        responsible: MT, EG and PW

class ExternalSolver 
 : public ImpuritySolver
 , public MatsubaraImpuritySolver 
{
public:
    /// @param executable the path to the executable
    ExternalSolver(const boost::filesystem::path& executable) 
    : exe_(boost::filesystem::complete(executable)) {}
    
    ImpuritySolver::result_type solve(
              const itime_green_function_t& G0
            , const alps::Parameters& parms);
    
    MatsubaraImpuritySolver::result_type solve_omega(
              const matsubara_green_function_t& G0_omega
            , const alps::Parameters& parms =alps::Parameters());
            
private:
    /// call the executable
    void call(std::string const& infile, std::string const& outfile);
      
    ///path to the solver executable
    boost::filesystem::path exe_;
};


#endif
