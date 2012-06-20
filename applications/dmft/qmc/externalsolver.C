/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>
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

/* $Id: externalsolver.C 360 2009-06-01 02:32:00Z gullc $ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver
#include "externalsolver.h"
#include "fouriertransform.h"
#include <cstdlib>
#include <fstream>
#ifdef BOOST_MSVC
#include <io.h>
#endif
#include "boost/tuple/tuple.hpp"
#include "alps/parser/parser.h"
#include "alps/utility/vectorio.hpp"
#include "alps/utility/temporary_filename.hpp"
#include "alps/numeric/isnan.hpp"
#include "alps/numeric/isinf.hpp"
void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta, const shape_t shape);
ImpuritySolver::result_type ExternalSolver::solve(const itime_green_function_t& G0, const alps::Parameters& parms) 
{
  alps::Parameters p(parms);
  //   BOOST_ASSERT(alps::is_master());
  std::string infile;
  std::string outfile; 
  if(p.defined("TMPNAME")){
    infile=p["TMPNAME"]+std::string(".in.h5");
    outfile=p["TMPNAME"]+std::string(".out.h5");
  } else{
    infile= alps::temporary_filename("alps_external_solver_in_");
    outfile= alps::temporary_filename("alps_external_solver_out_");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=infile;
  
  // write input file
  {
    alps::hdf5::archive solver_input(infile, "a");
    solver_input<<alps::make_pvp("/parameters",p);
    G0.write_hdf5(solver_input, "/G0");
  } 
  
  call(infile,outfile);
  
  // read the output
  unsigned int N=(unsigned int)p["N"];
  unsigned int sites     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int flavors   =(unsigned int)p.value_or_default("FLAVORS", 2);
  itime_green_function_t g(N+1, sites, flavors);
  {
    alps::hdf5::archive ar(outfile, "r");
    g.read_hdf5(ar, "/G0");
  }
  boost::filesystem::remove(boost::filesystem::absolute(outfile));
  return g;
}

MatsubaraImpuritySolver::result_type ExternalSolver::solve_omega(const matsubara_green_function_t& G0_omega, const alps::Parameters &parms)
{
  alps::Parameters p(parms);
  std::string infile;
  std::string outfile;
  if(p.defined("TMPNAME")){
    infile=p["TMPNAME"]+std::string(".in.h5");
    outfile=p["TMPNAME"]+std::string(".out.h5");
  } else{
    infile= alps::temporary_filename("alps_external_solver_in_");
    outfile= alps::temporary_filename("alps_external_solver_out_");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=outfile;
  {
    alps::hdf5::archive solver_input(infile, "a");
    if(p.value_or_default("SC_WRITE_DELTA", false)){
      //write Delta(i\omega_n) along with \Delta(\tau)
      
      //find the second moment of the band structure
      double epssqav ;
      if (parms.defined("DOSFILE")) {
        if (!parms.defined("EPSSQAV")) {
          throw std::logic_error("error: you specify a DOS file, please also specify the second moment of the band structure EPSSQAV!");
        } else {
          epssqav = parms["EPSSQAV"];
        }
      }else{
        double t=parms["t"]; //this is the default: semicircular density of states
        epssqav = t * t; //...and its moment.
      }
      double beta=p["BETA"];
      double n_orbital=p["FLAVORS"];
      double mu=p["MU"];
      double N=p["NMATSUBARA"];
      double N_tau=p["N"];
      double U=p["U"];
      FFunctionFourierTransformer Fourier(beta , 0, epssqav , n_orbital, 1);
      matsubara_green_function_t Delta_matsubara(N, 1, n_orbital);
      itime_green_function_t Delta_itime(N_tau+1, 1, n_orbital);
      for (int f = 0; f < n_orbital; ++f) {
        for (int i = 0; i < N; ++i) {
          Delta_matsubara(i, f) = -1. / G0_omega(i, f) + (std::complex < double >(mu, (2. * i + 1) * M_PI / beta));
        }
      }
      Fourier.backward_ft(Delta_itime, Delta_matsubara);
      Delta_matsubara.write_hdf5(solver_input, "/Delta_omega");
      Delta_itime.write_hdf5(solver_input, "/Delta");
      std::ofstream Delta_matsubara_file("Im_delta_matsubara.dat"); print_imag_green_matsubara(Delta_matsubara_file, Delta_matsubara, beta, diagonal);
      std::ofstream Delta_itime_file("Delta_itime.dat"); print_green_itime(Delta_itime_file,Delta_itime, beta, diagonal);
      
      //changed convention in the hybridization solvers:
      p["N_MATSUBARA"]=p["NMATSUBARA"];
      p["N_TAU"]=p["N"];
      p["N_ORBITALS"]=p["FLAVORS"];
      p["DELTA_IN_HDF5"]=1;
      p["DELTA"]=infile;
      p["MU"]=mu+U/2.;
    }
    solver_input<<alps::make_pvp("/parameters", p);
    G0_omega.write_hdf5(solver_input, "/G0");
  }
  
  call(infile,outfile);
  
  unsigned int n_matsubara=(unsigned int)p["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)p["N"];
  unsigned int n_site     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)p.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital);
  itime_green_function_t G_tau(n_tau+1, n_site, n_orbital);
  {
    alps::hdf5::archive ar(outfile, "r");
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  }
  //this is a safety check for impurity solvers.
  for(std::size_t i=0;i<n_orbital; ++i){
    for(std::size_t j=0;j<n_site;++j){
      for(std::size_t k=0;k<n_site;++k){
        for(std::size_t l=0;l<n_tau+1;++l){ 
          if(alps::numeric::isnan(G_tau(l,j,k,i)) || alps::numeric::isinf(G_tau(l,j,k,i))) {
            std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<i<<std::endl;
            std::cerr<<G_tau(l,j,k,i)<<" "<<std::endl;
            throw std::runtime_error("returned imag time Green's function contains nan or inf.");
          }
        }
        for(std::size_t l=0;l<n_matsubara;++l){ 
          if(alps::numeric::isnan(G_omega(l,j,k,i).real()) || alps::numeric::isnan(G_omega(l,j,k,i).imag())|| alps::numeric::isinf(G_omega(l,j,k,i).real()) || alps::numeric::isinf(G_omega(l,j,k,i).imag())) {
            std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<i<<std::endl;
            std::cerr<<G_omega(l,j,k,i).real()<<" "<<G_omega(l,j,k,i).imag()<<std::endl;
            throw std::runtime_error("returned freq Green's function contains nan or inf.");
          }
        }
      }
    }
  }
  boost::filesystem::remove(boost::filesystem::absolute(outfile));
  return std::make_pair(G_omega, G_tau);
}


void ExternalSolver::call(std::string const& infile, std::string const& outfile)
{
  
  // call the external solver program
  std::string command = exe_.string() + " " + infile + " " + outfile;
  //std::string command = "\""+exe_.string() + "\" " + infile + " " + outfile;
  //the line above won't work if the exe string has spaces in it.
  std::cerr << "Calling external solver " << exe_.string() << " as: "<<command<<std::endl;
  int result = std::system(command.c_str());
  if (result)
    boost::throw_exception(std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command));
  
  boost::filesystem::remove(infile);
  
  if (!boost::filesystem::exists(outfile))
    boost::throw_exception(std::runtime_error("The external impurity solver failed to write the output file named " + outfile));
}

