/*******************************************************************************
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#pragma once
#include <cstdlib>
#include <cstdio>
#include <string>
#include "alps/utilities/temporary_filename.hpp"
#include "alps/params.hpp"
#include "alps/gf/fourier.hpp"


class ImpuritySolver
{
public:
  ImpuritySolver(const std::string & launch_command, const std::string &param_file_name, const alps::params &p) :
    launch_command_(launch_command),
    param_file_name_(param_file_name),
    infile_str_(p["solver.INFILE_H5GF"].as<std::string>()),
    outfile_str_(p["solver.OUTFILE_H5GF"].as<std::string>()),
    params_(p),
    mu_(p["MU"]){}

  static void define_parameters(alps::params &p){
    p.define<std::string>("solver.INFILE_H5GF", "alps_solver_in.h5gf","Filename for temporary file");
    p.define<std::string>("solver.OUTFILE_H5GF", "alps_solver_out.h5gf","Filename for temporary file");
  }
  template<typename CTF> void write_cluster_info(const CTF& clustertransformer){
    alps::hdf5::archive solver_input(infile_str_, "w");
    clustertransformer.write_hdf5(solver_input, "/cluster");
  }
  template<typename MGF, typename IGF> void solve(const MGF& G0_omega, const IGF &G0_tau, MGF &G_omega, IGF &G_tau){
    std::cout<<"entering solve() command."<<std::endl;
    {
      alps::hdf5::archive solver_input(infile_str_, "w");
      solver_input<<alps::make_pvp("/parameters", params_);
      G0_omega.save(solver_input,"/G0");
      //also create an imaginary time version of G0
      G0_tau.save(solver_input,"/G0_tau");
    }
    call_solver();
    {
      alps::hdf5::archive ar(outfile_str_, "r");
      G_omega.load(ar, "/G_omega");
      G_tau  .load(ar, "/G_tau");
      validate_nan_inf(G_omega);
      validate_nan_inf(G_tau);
    }
    std::remove(infile_str_.c_str());
    std::remove(outfile_str_.c_str());
  }
  void adjust_mu(double mu){
    double mu_old=params_["MU"];
    params_["MU"]=mu;
    std::cout<<"adjusting solver mu from: "<<mu_old<<" to: "<<mu<<std::endl;
    mu_=mu;
  }

private:
  template<typename GF> void validate_nan_inf(const GF &G){
    //this is a safety check for impurity solvers.
    for(const typename GF::value_type *it=G.data().data(); it !=G.data().data()+G.data().num_elements(); it++)
      if(std::isnan(std::abs(*it)) || std::isinf(std::abs(*it))) throw std::runtime_error("encountered nan or inf in GF");
  }
  /// call the executable
  void call_solver(){
    std::stringstream  command; command << launch_command_ << " " << param_file_name_<< " --MU="<<mu_;
    std::cerr << "Calling external solver as: " <<  command.str()<<std::endl;
    int result = std::system(command.str().c_str());
    if (result)
      throw std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command.str());
    {
    std::ifstream check(outfile_str_.c_str());
    if(!check.good())
      throw std::runtime_error("The external impurity solver failed to write the output file named " + outfile_str_);
    }
  }

  ///command to launch the solver executable
  const std::string launch_command_;
  const std::string param_file_name_;
  const std::string infile_str_;
  const std::string outfile_str_;
  alps::params params_;
  double mu_;
};

inline void define_solver_parameters(alps::params &p){ 
  p.define<std::string>("solver.SOLVER", "path to the solver executable");
  ImpuritySolver::define_parameters(p);
}
