/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2026 by Emanuel Gull <gull@phys.columbia.edu>
 *               2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/


/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "interaction_expansion_choice.h"

#include <alps/parameter.h>
#include <alps/parameter_input.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/utilities/citations.hpp>
#include <alps/cli.hpp>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>
#include <filesystem>

#include "interaction_expansion/interaction_expansion.hpp"

namespace {
/// The external solvers that take the hybridization function Delta instead of G0.
bool reads_hybridization_function(const std::string& program) {
  return program == "hybridization" || program == "ct-hyb-segment";
}
}



/// @brief The DMFT main program
///
/// The program takes (at the moment) no command line options and reads the parameters 
/// of the simulation from the standard input. Parameters are read in the short ALPS text format where
/// each line contains one parameter definition in the format name=value .
///
/// At the moment the bare imaginary time green's function is assumed to be present in G_input_up and
/// G_input_down.
///
/// The following parameters have special meaning at this time: 
/// 
/// If SOLVER is Hirsch-Fye, a Hirsch-Fye solver is used, otherwise the value of SOLVER is assumed
/// to be the name of an external executable to be used as impurity solver.
/// @todo decent input routines needed!


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    if (alps::handle_cli_information(argc, argv, "dmft", [&] {
      std::cout << "Usage: " << argv[0] << " parameter_file\n";
    })) return 0;
    if (argc != 2) throw std::invalid_argument("Expected one parameter file");
    alps::cli_mpi_guard mpi(argc, argv);
    alps::Parameters parms = alps::read_parameter_file(argv[1]);
    parms["BASENAME"]=std::string(argv[1]);
    if (alps::cli_is_master()) {
      std::cout << "ALPS DMFT framework for the single site impurity problem.\n\n";
      alps::print_copyright(std::cout, "dmft");
    }
    // set working directory
    std::filesystem::path p(static_cast<std::string>(parms["BASENAME"]));
    if (!p.parent_path().empty())
      std::filesystem::current_path(p.parent_path());
    // all files are now relative to the input's directory
    const std::string input = p.filename().string();
    // results go to BASENAME.h5, which is the input itself for parm.h5
    if (p.extension() == ".h5")
      p.replace_extension();
    parms["BASENAME"]=p.filename().string();
    const std::string basename = static_cast<std::string>(parms["BASENAME"]);

    // ATOMIC made the interaction-expansion solver (the built-in one, or the
    // interaction program as an external solver) sample the atomic limit and
    // ignore the bath in every iteration, while the loop built the self-energy
    // and the next bath from the bath G0 (#347). Other solvers ignored it. It
    // is no longer supported, in either loop or with any solver.
    if (parms.defined("ATOMIC"))
      throw std::invalid_argument("dmft no longer supports ATOMIC. With the interaction-expansion "
                                  "solver (SOLVER = \"Interaction Expansion\" or the interaction "
                                  "program) it made the solver ignore the bath in every iteration "
                                  "(#347); other solvers ignored it. Remove ATOMIC from the "
                                  "parameter file.");

    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      // SYMMETRIZATION averages the two spins of every Green's function, which
      // undoes any field and leaves the impurity and lattice sides inconsistent.
      if (parms.value_or_default("SYMMETRIZATION",false)
          && (static_cast<double>(parms.value_or_default("H",0.))!=0.
              || static_cast<double>(parms.value_or_default("H_INIT",0.))!=0.))
        throw std::invalid_argument("SYMMETRIZATION = true requires H = 0 and H_INIT = 0: symmetrization "
                                    "averages the spins, so a field cannot act. For a magnetic field set "
                                    "SYMMETRIZATION = false.");
      if(!parms.value_or_default("OMEGA_LOOP",false)){
        //...imaginary time tau
        // Hirsch-Fye and CT-INT, built in or external, read G0 in Matsubara
        // frequencies only.
        const std::string solver = parms["SOLVER"];
        const std::string program = std::filesystem::path(solver).filename().string();
        if (solver=="Hirsch-Fye" || solver=="Interaction Expansion"
            || program=="hirschfye" || program=="interaction")
          throw std::invalid_argument("SOLVER = " + solver
                                      + " needs the Matsubara loop; set OMEGA_LOOP = 1");
        if (solver=="Hybridization")
          throw std::invalid_argument("The internal hybridization solver has been replaced by a standalone hybridzation solver.\nPlease use the \'hybridization\' program");
        // The imaginary-time loop builds the next hybridization function from
        // the Bethe-lattice condition Delta(tau) = t^2 G(tau) (ExternalSolver),
        // so it ran the Bethe lattice whatever TWODBS or DOSFILE said (#388,
        // case 3). L and tprime only act together with TWODBS. The Matsubara
        // loop has the general Hilbert transform.
        for (const char* lattice : {"TWODBS", "DOSFILE"})
          if (parms.defined(lattice))
            throw std::invalid_argument(std::string("The imaginary-time loop (OMEGA_LOOP = 0) supports only the "
                                        "Bethe lattice (semicircular DOS), and would ignore ") + lattice
                                        + "; set OMEGA_LOOP = 1 for " + lattice + ".");
        // Other solvers got G0(tau) and needed an imaginary-time Hilbert
        // transform, which exists only for the Matsubara loop: the run failed
        // after the first solver run (SemicircleHilbertTransformer::operator()).
        if (!reads_hybridization_function(program)
            && !static_cast<bool>(parms.value_or_default("SC_WRITE_DELTA", false)))
          throw std::invalid_argument("SOLVER = " + solver + " needs the Matsubara loop; set OMEGA_LOOP = 1. "
                                      "The imaginary-time loop (OMEGA_LOOP = 0) supports only solvers that read "
                                      "the hybridization function (hybridization, ct-hyb-segment, or "
                                      "SC_WRITE_DELTA = 1).");
        SemicircleHilbertTransformer transform(parms);
        std::cout<<"using external solver: "<<solver<<std::endl;
        parms["SC_WRITE_DELTA"]=1; //we need the hybridization function for this solver
        boost::shared_ptr<ImpuritySolver> solver_ptr(new ExternalSolver(solver));
        F_selfconsistency_loop(parms, *solver_ptr, transform);
      }
      else {
        //perform self consistency loop in Matsubara frequency omega
        if (parms.value_or_default("ANTIFERROMAGNET",false) && parms.value_or_default("SYMMETRIZATION",false)) { 
          std::cerr<<"ERROR: incompatible parameters: ANTIFERROMAGNET = true and SYMMETRIZATION = true"<<std::endl
                   <<"  for simulation in PM phase set: SYMMETRIZATION = true, ANTIFERROMAGNET = false."<<std::endl
                   <<"  for simulation in FM phase set: SYMMETRIZATION = false, ANTIFERROMAGNET = false."<<std::endl
                   <<"  for simulation in AFM phase set: SYMMETRIZATION = false, ANTIFERROMAGNET = true."<<std::endl;
          throw std::logic_error("Incompatible parameters: ANTIFERROMAGNET = true, SYMMETRIZATION = true.");
        }
          
        alps::scheduler::BasicFactory<InteractionExpansionSim,HubbardInteractionExpansionRun> interaction_expansion_factory_ss;
        alps::scheduler::BasicFactory<InteractionExpansionSim,MultiBandDensityHubbardInteractionExpansionRun> interaction_expansion_factory_mbd;
        // The solver keeps a reference to its factory, so every factory has to
        // live as long as solver_ptr.
        alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> hirschfye_factory;

        boost::shared_ptr<FrequencySpaceHilbertTransformer> transform_ptr;
        if (parms.value_or_default("SEMICIRCLE_HILBERT",false)) {
          transform_ptr.reset(new SemicircleFSHilbertTransformer(parms));
        } else {
          transform_ptr.reset(new GeneralFSHilbertTransformer(parms));
        }
        
        // Exactly one impurity solver is selected. The CT-INT run type is
        // chosen by select_interaction_expansion (FLAVORS/SITES read as ints);
        // every other SOLVER value is handed to the external/hybridization
        // path. No solver is constructed-then-discarded.
        boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;
        if (parms["SOLVER"]=="Interaction Expansion") {
          const int flavors = static_cast<int>(parms.value_or_default("FLAVORS", 2));
          const int sites   = static_cast<int>(parms.value_or_default("SITES", 1));
          switch (select_interaction_expansion(flavors, sites)) {
            case interaction_expansion_choice::single_site_hubbard:
              std::cout<<"using single site Hubbard solver"<<std::endl;
              solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_ss,parms,argc,argv));
              break;
            case interaction_expansion_choice::multiband_density:
              std::cout<<"using multiband Hubbard solver"<<std::endl;
              solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_mbd,parms,argc,argv));
              break;
            case interaction_expansion_choice::unsupported:
              throw std::runtime_error("DMFT Interaction Expansion: unsupported (FLAVORS, SITES) "
                                       "combination; set FLAVORS=2 (single site) or SITES=1 (multiband).");
          }
          if (alps::cli_is_master()) alps::print_citations(std::cout, "interaction");
        }
        else if (parms["SOLVER"]=="Hirsch-Fye") {
          std::cout<<"solving Hirsch Fye"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(hirschfye_factory,parms,argc,argv));
          // Built-in solvers print no notice of their own; external solvers do.
          if (alps::cli_is_master()) alps::print_citations(std::cout, "hirschfye");
        }
        else if (parms["SOLVER"]=="Hybridization") {
          throw std::invalid_argument("The internal hybridization solver has been replaced by a standalone hybridzation solver.\nPlease use the \'hybridization\' program");
        }
        else {
          std::string p(parms["SOLVER"]);
          std::cout<<"using external solver: "<<p<<std::endl;
          if(reads_hybridization_function(std::filesystem::path(p).filename().string()))
            parms["SC_WRITE_DELTA"]=1; //we need the hybridization function for this solver
          solver_ptr.reset(new ExternalSolver(p));
        }
        selfconsistency_loop_omega(parms, *solver_ptr, *transform_ptr);
      }
    }
    else { //CLUSTER_LOOP
      throw std::logic_error("you should use the cluster framework for a cluster calculation.");
    }    
    // Rewrite the input parameter file in its own format. If it cannot be
    // rewritten (for example a read-only input file) warn instead of failing:
    // the results are intact and BASENAME.h5 is still written below.
    try {
      alps::write_parameter_file(input, parms);
    } catch (std::exception const& exc) {
      std::cerr << "Warning: could not rewrite parameter file '" << input
                << "': " << exc.what() << "\n";
    }
    if (basename+".h5" != input)
      alps::write_parameter_file(basename+".h5", parms);
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr<<exc.what()<<std::endl;
    return -1;
  }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif  
  return 0;  
}
