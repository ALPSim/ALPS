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

/* $Id: externalsolver.C 360 2009-06-01 02:32:00Z gullc $ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver
#include "externalsolver.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "bandstructure.h"
#include <cstdlib>
#include <fstream>
#ifdef BOOST_MSVC
#include <io.h>exe
#endif
#include "boost/tuple/tuple.hpp"
#include "alps/parser/parser.h"
#include "alps/utility/vectorio.hpp"
#include "alps/utility/temporary_filename.hpp"
#include "alps/numeric/isnan.hpp"
#include "alps/numeric/isinf.hpp"
#include <alps/utility/os.hpp>

ExternalSolver::ExternalSolver(const boost::filesystem::path& executable)
: exe_(boost::filesystem::absolute(executable,alps::bin_directory()))
{
}


void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta, const shape_t shape);

void prepare_parms_for_hybridization(alps::Parameters& p, alps::hdf5::archive& solver_input) {
  double mu=p["MU"];
  U_matrix U(p);
  double h = static_cast<double>(p.value_or_default("H",0.));
  unsigned n_orbital=static_cast<unsigned>(p["FLAVORS"]);
  
  //changed convention in the hybridization solvers:
  p["N_MATSUBARA"]=p["NMATSUBARA"];
  p["N_TAU"]=p["N"];
  p["N_ORBITALS"]=p["FLAVORS"];
  p["DELTA_IN_HDF5"]=1;
  p["DELTA"]=p["INFILE"];
  if (h==0) {
    p["MU"]=mu+U.mu_shift();
  } else {
    std::vector<double> mu_vector(n_orbital);
    for (int f = 0; f < n_orbital; ++f) {
      double hsign = f%2 ? h : -h;
      mu_vector[f]=mu+U.mu_shift()+hsign;
    }
    p["MU_VECTOR"]=p["INFILE"];
    p["MU_IN_HDF5"]=1;
    solver_input<<alps::make_pvp("/MUvector",mu_vector);
    p.erase("MU");
  }
  p["DMFT_FRAMEWORK"]=1;
}

ImpuritySolver::result_type ExternalSolver::solve(const itime_green_function_t& G0, const alps::Parameters& parms) 
{
  alps::Parameters p(parms);
  unsigned int n_tau=(unsigned int)parms["N"];
  unsigned int n_site     =(unsigned int)parms.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)parms.value_or_default("FLAVORS", 2);
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
    if(parms.value_or_default("SC_WRITE_DELTA", false)){
      //write \Delta(\tau)
      //note: G0 is the full G
      
      SemicircleBandstructure bethe_parm(parms);  // to get the second moment(s) of the band structure
      itime_green_function_t Delta_itime(n_tau+1, 1, n_orbital);
      bool AFM=parms.value_or_default("ANTIFERROMAGNET",false);
      for (int f = 0; f < n_orbital; ++f) {
        int f_=AFM?(f%2==0?f+1:f-1):f;
        for (int i = 0; i <= n_tau; ++i) {
          Delta_itime(i, f) = bethe_parm.second_moment(f) * G0(i, f_);
        }
      }
      Delta_itime.write_hdf5(solver_input, "/Delta");
      std::ofstream Delta_itime_file("Delta_itime.dat"); print_green_itime(Delta_itime_file,Delta_itime, static_cast<double>(p["BETA"]), diagonal);
      
      prepare_parms_for_hybridization(p,solver_input);
    } else {
      G0.write_hdf5(solver_input, "/G0");
    }
    solver_input<<alps::make_pvp("/parameters",p);
  } 
  
  call(infile,outfile);
  
  // read the output
  itime_green_function_t g(n_tau+1, n_site, n_orbital);
  {
    alps::hdf5::archive ar(outfile, "r");
    g.read_hdf5(ar, "/G_tau");
  }
  boost::filesystem::remove(outfile);
  return g;
}

MatsubaraImpuritySolver::result_type ExternalSolver::solve_omega(const matsubara_green_function_t& G0_omega, const alps::Parameters &parms)
{
  alps::Parameters p(parms);
  unsigned int n_matsubara=(unsigned int)parms["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)parms["N"];
  unsigned int n_site     =(unsigned int)parms.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)parms.value_or_default("FLAVORS", 2);
  std::string infile;
  std::string outfile;
  if(parms.defined("TMPNAME")){
    infile=parms["TMPNAME"]+std::string(".in.h5");
    outfile=parms["TMPNAME"]+std::string(".out.h5");
  } else{
    infile= alps::temporary_filename("alps_external_solver_in_");
    outfile= alps::temporary_filename("alps_external_solver_out_");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=outfile;
  {
    alps::hdf5::archive solver_input(infile, "a");
    if(parms.value_or_default("SC_WRITE_DELTA", false)){
      //write Delta(i\omega_n) along with \Delta(\tau)
      
      double beta=p["BETA"];
      double mu=p["MU"];
      double h = static_cast<double>(parms.value_or_default("H",0.));
      FFunctionFourierTransformer Fourier(parms);
      matsubara_green_function_t Delta_matsubara(n_matsubara, 1, n_orbital);
      itime_green_function_t Delta_itime(n_tau+1, 1, n_orbital);
      for (int f = 0; f < n_orbital; ++f) {
        double hsign = f%2 ? h : -h;
        for (int i = 0; i < n_matsubara; ++i) {
          Delta_matsubara(i, f) = -1. / G0_omega(i, f) + (std::complex < double >(mu+hsign, (2. * i + 1) * M_PI / beta));
        }
      }
      Fourier.backward_ft(Delta_itime, Delta_matsubara);
      Delta_matsubara.write_hdf5(solver_input, "/Delta_omega");
      Delta_itime.write_hdf5(solver_input, "/Delta");
      std::ofstream Delta_matsubara_file("Im_delta_matsubara.dat"); print_imag_green_matsubara(Delta_matsubara_file, Delta_matsubara, beta, diagonal);
      std::ofstream Delta_itime_file("Delta_itime.dat"); print_green_itime(Delta_itime_file,Delta_itime, beta, diagonal);
      
      prepare_parms_for_hybridization(p,solver_input);
    }
    solver_input<<alps::make_pvp("/parameters", p);//hier problem
    G0_omega.write_hdf5(solver_input, "/G0");
  }
  
  call(infile,outfile);
  
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital);
  itime_green_function_t G_tau(n_tau+1, n_site, n_orbital);
  alps::hdf5::archive ar(outfile, "r");
  if (ar.is_group("G_omega")) {    
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  } else { // Fourier: G_tau -> G_omega
    G_tau.read_hdf5(ar, "/G_tau");
    boost::shared_ptr < FourierTransformer > fourier_ptr;
    std::cout<<"converting ALPS parameters"<<std::endl;
    std::vector<double>n(n_orbital,0.);
    for(unsigned int u=0;u<n_orbital;u++){
      n[u]=-G_tau(n_tau,0,0,u);  // convention used in the selfconsistency: G(tau=beta^-,orbital)=-density(orbital); same on the output of hybridization solver
    }
    FourierTransformer::generate_transformer_U(parms, fourier_ptr, n); //still takes old alps parameter class.
    fourier_ptr->forward_ft(G_tau,G_omega);
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
  boost::filesystem::remove(outfile);
  return std::make_pair(G_omega, G_tau);
}


void ExternalSolver::call(std::string const& infile, std::string const& outfile)
{
  
  // call the external solver program
  std::string command = "\""+exe_.string() + "\" " + infile + " " + outfile; 
  //the line above won't work if the exe string has spaces in it.
  //we need to use quotes if it has spaces - where does it not work?
  std::cerr << "Calling external solver " << exe_.string() << " as: "<<command<<std::endl;
  int result = std::system(command.c_str());
  if (result)
    boost::throw_exception(std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command));
  
  boost::filesystem::remove(infile);
  
  if (!boost::filesystem::exists(outfile))
    boost::throw_exception(std::runtime_error("The external impurity solver failed to write the output file named " + outfile));
}

