/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "maxent.hpp"

/*

alps::scheduler::Task* MaxEntFactory::make_task(const alps::ProcessList& w, 
                 const boost::filesystem::path& fn) const
{
  return static_cast<alps::scheduler::Task*>(new MaxEntSimulation(w,fn));
}
  



void MaxEntFactory::print_copyright(std::ostream& out) const
{
  out << "ALPS Maximum Entropy application\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>\n"
      << "                        Thomas Pruschke <pruschke@comp-phys.org>\n"
      << "                        Matthias Troyer <troyer@comp-phys.org>\n"
      << "  copyright (c) 2012 by Emanuel Gull <gull@phys.columbia.edu>\n"
      << " for details see the publication:\n"
      << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}

*/

bool stop_callback(boost::posix_time::ptime const & end_time) {
    static alps::ngs::signal signal;
    return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}


#ifdef BUILD_PYTHON_MODULE
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
namespace nb = nanobind;

void run_it(nb::dict const & parms_){
  // MaxEntSimulation derives from alps::mcbase, and stop_callback above
  // constructs an alps::ngs::signal, so running it installs ALPS handlers
  // for SIGINT and for SIGSEGV/SIGBUS. Restore Python's before returning
  // to the interpreter, exactly as cthyb::solve and ctint::solve do.
  pyalps::scoped_signal_handlers signal_handlers;
    alps::parameters_type<MaxEntSimulation>::type parms = pyalps::params_from_dict(parms_);
    std::string out_file = boost::lexical_cast<std::string>(parms["BASENAME"]|"results")+std::string(".out.h5");

#else

  int main(int argc, char** argv)
{
  alps::mcoptions options(argc, argv);
    
  alps::parameters_type<MaxEntSimulation>::type parms(alps::hdf5::archive(options.input_file));
  
    std::string out_file(boost::lexical_cast<std::string>(parms["BASENAME"]|options.output_file)+std::string(".out.h5"));
//    std::cout << out_file << std::endl;
#endif
  MaxEntSimulation my_sim(parms,out_file); // creat a simulation
  my_sim.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)(parms["MAX_TIME"]|60)))); // run the simulation
#ifdef BUILD_PYTHON_MODULE
  return;
#else
/*
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    return alps::scheduler::start(argc, argv, MaxEntFactory());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
    return -1;
  }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
*/
  return 0;
#endif
}
    
#ifdef BUILD_PYTHON_MODULE
    NB_MODULE(maxent_c, m) {
        m.def("AnalyticContinuation", run_it);
    }
#endif
