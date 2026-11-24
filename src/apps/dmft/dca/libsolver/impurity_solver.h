/*******************************************************************************
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#pragma once
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "alps/hdf5.hpp"
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
  template<typename MGF, typename IGF> void solve(const MGF& G0_omega, const IGF &G0_tau, MGF &G_omega, IGF &G_tau){
    std::cout<<"entering solve() command."<<std::endl;
    // Only a file the solver writes now can be read as its result: an output
    // file left from before would be taken for it if the solver exits
    // normally without writing one (#459). It is removed before the input is
    // written, which may be the same file.
    std::remove(outfile_str_.c_str());
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
      // load() takes the meshes from the file. G_tau is only printed, on
      // whatever tau grid the solver used.
      check_meshes(G_omega, G0_omega);
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
  /// The solver's G_omega has to be on the meshes of the G0 it was given:
  /// another BETA, NMATSUBARA, cluster or number of flavors would be combined
  /// with G0 on the loop's meshes (#459).
  template<typename GF> void check_meshes(const GF &G, const GF &G0) const{
    typedef typename GF::gf_type gf_type;
    if (!(static_cast<const gf_type&>(G).meshes() == static_cast<const gf_type&>(G0).meshes()))
      throw std::runtime_error(outfile_str_ + ": the solver's G_omega is not on the meshes of G0"
                               " (BETA, NMATSUBARA, cluster, FLAVORS) that dca gave it.");
  }
  template<typename GF> void validate_nan_inf(const GF &G){
    //this is a safety check for impurity solvers.
    for(const typename GF::value_type *it=G.data().data(); it !=G.data().data()+G.data().num_elements(); it++)
      if(std::isnan(std::abs(*it)) || std::isinf(std::abs(*it))) throw std::runtime_error("encountered nan or inf in GF");
  }
  /// x with the fewest digits that read back as x (0.3, not 0.29999999999999999)
  static std::string exact(double x){
    for (int digits = 1; digits < std::numeric_limits<double>::max_digits10; ++digits) {
      std::ostringstream out;
      out.precision(digits);
      out << x;
      if (std::stod(out.str()) == x) return out.str();
    }
    std::ostringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << x;
    return out.str();
  }
  /// s as one word for the shell: in single quotes, with each ' as '\''
  static std::string shell_quote(const std::string &s){
    std::string quoted = "'";
    for (char c : s) quoted += (c == '\'' ? std::string("'\\''") : std::string(1, c));
    return quoted + "'";
  }
  /// call the executable
  void call_solver(){
    // solver.SOLVER is a command line (it may hold a launcher such as mpirun);
    // the parameter file is one argument, quoted for the shell, and MU has the
    // digits that read back as the same double (#460).
    std::stringstream  command;
    command << launch_command_ << " " << shell_quote(param_file_name_) << " --MU=" << exact(mu_);
    std::cerr << "Calling external solver as: " <<  command.str()<<std::endl;
    int result = std::system(command.str().c_str());
    if (result)
      throw std::runtime_error("System error code " +std::to_string(result) + " encountered when executing command:\n"+command.str());
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
