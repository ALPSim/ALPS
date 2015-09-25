/*******************************************************************************
 *
 * Copyright (C) 2005-2015  Emanuel Gull <gullc@comp-phys.org>
 *
 *****************************************************************************/

#pragma once
#include <cstdlib>
#include "alps/utilities/temporary_filename.hpp"

class ImpuritySolver
{
public:
  ImpuritySolver(const std::string & launch_command, const alps::params &p) :
    launch_command_(launch_command),
    params_(p),
    infile_str_(p["INFILE"].as<std::string>()),
    outfile__str_(p["INFILE"]){}

  static void define_parameters(alps::params &p){
    p.define<std::string>("INFILE", alps::temporary_filename("alps_external_solver_")+".in.h5","Filename for temporary file");
    p.define<std::string>("OUTFILE", alps::temporary_filename("alps_external_solver_")+".out.h5","Filename for temporary file");
  }
  template<typename CTF> void write_cluster_info(const CTF& clustertransformer){
    alps::hdf5::archive solver_input(infile_str_, alps::hdf5::archive::WRITE);
    clustertransformer.write_hdf5(solver_input, "/cluster");
  }
  template<typename MGF, typename IGF> void solve(const MGF& G0_omega, MGF &G_omega, IGF &G_tau){
    {
      alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
      solver_input<<alps::make_pvp("/parameters", params_);
      solver_input<<alps::make_pvp("/G0",G0_omega);
    }
    call_solver();


    {
      alps::hdf5::archive ar(outfile_str_, alps::hdf5::archive::READ);
      ar>>alps::make_pvp(G_omega, "/G_omega");
      ar>>alps::make_pvp(G_tau  , "/G_tau");
      validate_nan_inf(G_omega);
      validate_nan_inf(G_tau);
    }
    boost::filesystem::remove(infile_str_);
    boost::filesystem::remove(outfile_str_);
  }

private:
  template<typename GF> void validate_nan_inf(const GF &G){
    //this is a safety check for impurity solvers.
    for(typename GF::container_type::const_iterator it=G.data().begin(); it!=G.data().end();++it){
      if(std::isnan(*it) || std::isinf(*it) ){
        std::cerr<<G<<std::endl;
        throw std::runtime_error("Green's function obtained from solver contains nan or inf.");
      }
    }
  }
  /// call the executable
  void call_solver(){

    // call the external solver program
    std::string command = launch_command_ + " " + infile_str_ + " " + outfile_str_;
    std::cerr << "Calling external solver as: " <<  command<<std::endl;
    int result = std::system(command.c_str());
    if (result)
      throw std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command);

    if (!boost::filesystem::exists(outfile_str_))
      throw std::runtime_error("The external impurity solver failed to write the output file named " + outfile_str_);
  }


  ///command to launch the solver executable
  const std::string launch_command_;
  const std::string infile_str_;
  const std::string outfile_str_;
  const alps::params params_;
};
