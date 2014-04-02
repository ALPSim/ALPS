#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
#include <unistd.h>
#include <boost/program_options.hpp>


int main(int argc, char **argv){
  namespace po = boost::program_options;
  // Declare the supported options.
  po::options_description desc("Allowed options");
  std::string job_type;
  std::string job_class;
  std::string executable;
  std::string bg_size;
  std::string initialdir;
  std::string input;
  std::string output;
  std::string error;
  std::string wall_clock_limit;
  std::string notification;
  std::string notify_user;
  std::stringstream arguments;
  std::string exe;
  std::string mode;
  std::string args;
  std::string infile_name;
  std::string outfile_name;
  std::string jobfile_name;
 
  char buf[1024]; 
desc.add_options()
    ("help", "show this help")
    ("job_type", po::value<std::string>(&job_type)->default_value("bluegene"), "job type value, should be bluegene")
    ("job_class", po::value<std::string>(&job_class)->default_value("short"), "short, normaldyn, or long.")
    ("executable", po::value<std::string>(&executable)->default_value("/bgl/BlueLight/ppcfloor/bglsys/bin/mpirun"), "executable: this should be the mpirun program")
    ("bg_size", po::value<std::string>(&bg_size)->default_value("1024"), "32, 128 (for long), 512 (for normaldyn), and 512, 1024, ... for short")
    ("initialdir", po::value<std::string>(&initialdir)->default_value(getcwd(buf, 1024)), "job run directory (current directory detected)")
    ("input", po::value<std::string>(&input)->default_value("/dev/null"), "program std input")
    ("output", po::value<std::string>(&output)->default_value("$(jobid).out"), "program std output")
    ("output", po::value<std::string>(&error)->default_value("$(jobid).err"), "program std error")
    ("wall_clock_limit", po::value<std::string>(&wall_clock_limit)->default_value("01:00:00"), "program wall clock time, format hh:mm:ss")
    ("notification", po::value<std::string>(&notification)->default_value("complete"), "do we notify user after completion?")
    ("notify_user", po::value<std::string>(&notify_user)->default_value("gull@phys.columbia.edu"), "e-mail address of user")
    ("exe", po::value<std::string>(&exe)->default_value("/gpfs/home3/g/gull/alps_svn/DMFT-solvers/CTAUX_SUB_BG/ctaux_mpi"), "actual executable run by mpirun")
    ("mode", po::value<std::string>(&mode)->default_value("CO"), "coprocessor CO or virtual node VN mode")
    ("infile_name", po::value<std::string>(&infile_name), "input filename")
    ("outfile_name", po::value<std::string>(&outfile_name), "output filename")
    //("jobfile_name", po::value<std::string>(&jobfile_name)->default_value("ctaux.run"), "job filename")
          ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }
  if(!vm.count("infile_name")) throw std::runtime_error("you need to specify the parameter file as infile_name");
  if(!vm.count("outfile_name")) throw std::runtime_error("you need to specify the parameter file as outfile_name");

  std::string prefix("# @ ");

  std::cout<<prefix<<"job_type"<<" = "<<job_type<<std::endl;
  std::cout<<prefix<<"class"<<" = "<<job_class<<std::endl;
  std::cout<<prefix<<"executable"<<" = "<<executable<<std::endl;
  std::cout<<prefix<<"bg_size"<<" = "<<bg_size<<std::endl;
  std::cout<<prefix<<"initialdir"<<" = "<<initialdir<<std::endl;
  std::cout<<prefix<<"input"<<" = "<<input<<std::endl;
  std::cout<<prefix<<"output"<<" = "<<output<<std::endl;
  std::cout<<prefix<<"error"<<" = "<<error<<std::endl;
  std::cout<<prefix<<"wall_clock_limit"<<" = "<<wall_clock_limit<<std::endl;
  std::cout<<prefix<<"notification"<<" = "<<notification<<std::endl;
  std::cout<<prefix<<"notify_user"<<" = "<<notify_user<<std::endl;

  arguments<<prefix<<"arguments"<<" = "<<" -exe "<<exe<<" -mode "<<mode<<" -args "<<"\""<<infile_name<<" "<<outfile_name<<"\""<<" -cwd "<<initialdir;
  std::cout<<arguments.str()<<std::endl;
  std::cout<<prefix<<" queue "<<std::endl;
}

