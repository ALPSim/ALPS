/*******************************************************************************
 *
 * Copyright (C) 2005-2015  Emanuel Gull <gullc@comp-phys.org>
 *
 *****************************************************************************/

#pragma once
#include <cstdlib>
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
    params_(p){}

  static void define_parameters(alps::params &p){
    p.define<std::string>("solver.INFILE_H5GF", "alps_solver_in.h5gf","Filename for temporary file");
    p.define<std::string>("solver.OUTFILE_H5GF", "alps_solver_out.h5gf","Filename for temporary file");
  }
  template<typename CTF> void write_cluster_info(const CTF& clustertransformer){
    alps::hdf5::archive solver_input(infile_str_, alps::hdf5::archive::WRITE);
    clustertransformer.write_hdf5(solver_input, "/cluster");
  }
  template<typename MGF, typename IGF> void solve(const MGF& G0_omega, const IGF &G0_tau, MGF &G_omega, IGF &G_tau){
    {
      alps::hdf5::archive solver_input(infile_str_, alps::hdf5::archive::WRITE);
      solver_input<<alps::make_pvp("/parameters", params_);
      G0_omega.save(solver_input,"/G0");
      //also create an imaginary time version of G0
      G0_tau.save(solver_input,"/G0_tau");
    }
    call_solver();
    {
      alps::hdf5::archive ar(outfile_str_, alps::hdf5::archive::READ);
      G_omega.load(ar, "/G_omega");
      G_tau  .load(ar, "/G_tau");
      validate_nan_inf(G_omega);
      validate_nan_inf(G_tau);
    }
    boost::filesystem::remove(infile_str_);
    boost::filesystem::remove(outfile_str_);
  }

private:
  template<typename GF> void validate_nan_inf(const GF &G){
    //this is a safety check for impurity solvers.
    for(const typename GF::value_type *it=G.data().data(); it !=G.data().data()+G.data().num_elements(); it++)
      if(std::isnan(std::abs(*it)) || std::isinf(std::abs(*it))) throw std::runtime_error("encountered nan or inf in GF");
  }
  /// call the executable
  void call_solver(){
    std::string command = launch_command_ + " " + param_file_name_;
    std::cerr << "Calling external solver as: " <<  command<<std::endl;
    int result = std::system(command.c_str());
    if (result)
      throw std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command);

    if (!boost::filesystem::exists(outfile_str_))
      throw std::runtime_error("The external impurity solver failed to write the output file named " + outfile_str_);
  }

  ///command to launch the solver executable
  const std::string launch_command_;
  const std::string param_file_name_;
  const std::string infile_str_;
  const std::string outfile_str_;
  const alps::params params_;
};

inline void define_solver_parameters(alps::params &p){ 
  p.define<std::string>("solver.SOLVER", "path to the solver executable");
  ImpuritySolver::define_parameters(p);
}
