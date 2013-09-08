/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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
//compile it as a python module (requires boost::python library)
using namespace boost::python;

void run_it(boost::python::dict parms_){
    alps::parameters_type<MaxEntSimulation>::type parms(parms_);
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
    BOOST_PYTHON_MODULE(maxent_c)
    {
        def("AnalyticContinuation",run_it);//define python-callable run method
    };
#endif

