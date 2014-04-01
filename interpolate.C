/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2011 by Emanuel Gull <gull@phys.columbia.edu>
 * Copyright (C) 2009 by Sebastian Fuchs <fuchs@comp-phys.org>
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

#include "interpolator.h"
#include "cluster.h"
#include <alps/alea.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;



int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    po::options_description desc("Allowed options");
    std::string pfilename;
    std::string selfenergy_file;
    desc.add_options()
    ("help", "show this help")
    ("pfile", po::value<std::string>(&pfilename)->default_value("nambu.param"), "parameter file")
    ("selfenergy_file", po::value<std::string>(&selfenergy_file)->default_value(""), "selfenergy file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    
    if (vm.count("help")) {
      std::cout<<desc;
      return 1;
    }
    std::ifstream parm_file(pfilename.c_str());
    if(!parm_file.is_open()) throw std::runtime_error("parm file not open.");
    alps::Parameters parms(parm_file);


    ClusterTransformer* clusterhandler=ClusterTransformer::generate_transformer(parms);
    
    std::vector<double> values(8,1.);
    std::vector<double> kvec(2); kvec[0]=0.5; kvec[1]=0.5;
    double z=clusterhandler->interpolate(kvec,values, false);
    for(int i=0;i<8;++i) values[i]=i;
    
    int N=10;
    for(int i=-N;i<N;++i){
      for(int j=-N;j<N;++j){
        double kx=M_PI/N*i;
        double ky=M_PI/N*j;
        kvec[0]=kx; kvec[1]=ky;
        std::cout<<kx<<" "<<ky<<" "<<clusterhandler->interpolate(kvec,values, false)<<std::endl;
      }
    }
    delete clusterhandler;
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << "\n";
    std::exit(-5);
  }
#endif
}



