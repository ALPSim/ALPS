#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
#include <boost/program_options.hpp>
#include <alps/hdf5.hpp>


int main(int argc, char **argv){
  namespace po = boost::program_options;
  std::string new_K;
  std::string new_KMAX;
  std::string new_MAX_TIME;
  std::string infile;
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show this help")
    ("K", po::value<std::string>(&new_K), "change the value of the ct-aux parameter K")
    ("KMAX", po::value<std::string>(&new_KMAX), "change the value of the ct-aux parameter KMAX")
    ("MAX_TIME", po::value<std::string>(&new_MAX_TIME), "change the value of MAX_TIME")
    ("h5", po::value<std::string>(&infile), "file to be changed")
          ;
  po::variables_map vm;
  po::positional_options_description p;
  p.add("h5", 1);

  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }

  alps::hdf5::archive arch(infile, alps::hdf5::archive::WRITE);
  if(new_K != ""){
    std::cout<<"resetting K to: "<<new_K<<std::endl;
    arch<<alps::make_pvp("/parameters/K",new_K);
    std::cout<<"done resetting "<<std::endl;
  }
  if(new_KMAX != ""){
    std::cout<<"resetting KMAX to: "<<new_KMAX<<std::endl;
    arch<<alps::make_pvp("/parameters/KMAX",new_KMAX);
  }
  if(new_MAX_TIME != ""){
    std::cout<<"resetting new MAX_TIME to: "<<new_MAX_TIME<<std::endl;
    arch<<alps::make_pvp("/parameters/KMAX",new_MAX_TIME);
  }
}

