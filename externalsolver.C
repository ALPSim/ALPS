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

/* $Id$ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver

#include <cstdlib>
#include <fstream>

#include "alps/utility/temporary_filename.hpp"
#include "externalsolver.h"
#include "alps/parameter.h"
#include "cluster.h"

#include <boost/lexical_cast.hpp>

ImpuritySolver::matsubara_result_type ExternalSolver::solve_matsubara(
                                                                 const matsubara_green_function_t& G0_omega
                                                                 , alps::Parameters p, ClusterTransformer *c_ptr)
{
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
    alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
    solver_input<<alps::make_pvp("/parameters", p);
    G0_omega.write_hdf5(solver_input, "/G0");
    
    if((bool)(p.value_or_default("WRITE_CLUSTER_INFO",false))){
      if(c_ptr==0) throw std::runtime_error("called cluster info function with a solver equal to 0");
      c_ptr->write_hdf5(solver_input, "/cluster");
    }
  }
  
  call(infile,outfile);
  
  
  unsigned int n_matsubara=(unsigned int)p["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)p["N"];
  unsigned int n_site     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)p.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital, G0_omega.shape());
  itime_green_function_t G_tau(n_tau+1, n_site, n_orbital, G0_omega.shape());
  {
    alps::hdf5::archive ar(outfile, alps::hdf5::archive::READ);
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  }
  //this is a safety check for impurity solvers.
  for(std::size_t i=0;i<n_orbital; ++i){
    for(std::size_t j=0;j<n_site;++j){
      for(std::size_t k=0;k<n_site;++k){
        for(std::size_t l=0;l<n_tau+1;++l){ 
          if(std::isnan(G_tau(l,j,k,i)) || std::isinf(G_tau(l,j,k,i))) {
            std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<i<<std::endl;
            std::cerr<<G_tau(l,j,k,i)<<" "<<std::endl;
            throw std::runtime_error("returned imag time Green's function contains nan or inf.");
          }
        }
        for(std::size_t l=0;l<n_matsubara;++l){ 
          if(std::isnan(G_omega(l,j,k,i).real()) || std::isnan(G_omega(l,j,k,i).imag())|| std::isinf(G_omega(l,j,k,i).real()) || std::isinf(G_omega(l,j,k,i).imag())) {
            std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<i<<std::endl;
            std::cerr<<G_omega(l,j,k,i).real()<<" "<<G_omega(l,j,k,i).imag()<<std::endl;
            throw std::runtime_error("returned freq Green's function contains nan or inf.");
          }
        }
      }
    }
  }
  boost::filesystem::remove(boost::filesystem::system_complete(outfile));
  return std::make_pair(G_omega, G_tau);
}

ImpuritySolver::nambu_result_type ExternalSolver::solve_nambu(
                                                             const nambu_matsubara_green_function_t& G0_omega
                                                             , alps::Parameters p)
{
  std::cout<<"solve nambu for external solver!"<<std::endl;
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
    alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
    solver_input<<alps::make_pvp("/parameters", p);
    G0_omega.write_hdf5(solver_input, "/G0");
  }
  
  call(infile,outfile);
  
  unsigned int n_matsubara=(unsigned int)p["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)p["N"];
  unsigned int n_site     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)p.value_or_default("FLAVORS", 2);
  nambu_matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital, nondiagonal);
  nambu_itime_green_function_t G_tau(n_tau+1, n_site, n_orbital, nondiagonal);
  {
    alps::hdf5::archive ar(outfile, alps::hdf5::archive::READ);
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  }
  //this is a safety check for impurity solvers.
  for(std::size_t f1=0;f1<n_orbital; ++f1){
    for(std::size_t f2=0;f2<n_orbital; ++f2){
      for(std::size_t j=0;j<n_site;++j){
        for(std::size_t k=0;k<n_site;++k){
          if(j!=k && (G_tau.shape()==diagonal)) continue;
          for(std::size_t l=0;l<n_tau+1;++l){
            if(std::isnan(G_tau(l,j,k,f1,f2)) || std::isinf(G_tau(l,j,k,f1,f2))) {
              std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<f1<<" "<<f2<<std::endl;
              std::cerr<<G_tau(l,j,k,f1,f2)<<" "<<std::endl;
              throw std::runtime_error("returned imag time Green's function contains nan or inf.");
            }
          }
          for(std::size_t l=0;l<n_matsubara;++l){ 
            if(std::isnan(G_omega(l,j,k,f1,f2).real()) || std::isnan(G_omega(l,j,k,f1,f2).imag())|| std::isinf(G_omega(l,j,k,f1,f2).real()) || std::isinf(G_omega(l,j,k,f1,f2).imag())) {
              std::cerr<<"freq: "<<l<<" sites: "<<j<<" "<<k<<" spin: "<<f1<<" "<<f2<<std::endl;
              std::cerr<<G_omega(l,j,k,f1,f2).real()<<" "<<G_omega(l,j,k,f1,f2).imag()<<std::endl;
              throw std::runtime_error("returned freq Green's function contains nan or inf.");
            }
          }
        }
      }
    }
  }
  boost::filesystem::remove(boost::filesystem::system_complete(outfile));
  return std::make_pair(G_omega, G_tau);
}


void ExternalSolver::call(std::string const& infile, std::string const& outfile)
{
  
  // call the external solver program
  std::string command = exe_.native_file_string() + " " + infile + " " + outfile;
  std::cerr << "Calling external solver " << exe_.native_file_string() << " as: "<<command<<std::endl;
  int result = std::system(command.c_str());
  if (result)
    throw std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command);
  
  boost::filesystem::remove(infile);
  
  if (!boost::filesystem::exists(outfile))
    throw std::runtime_error("The external impurity solver failed to write the output file named " + outfile);
}
