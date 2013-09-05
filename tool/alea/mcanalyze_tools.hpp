/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef CPP_COMMANDLINE_HPP
#define CPP_COMMANDLINE_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/utility/encode.hpp>
#include <alps/numeric/vector_functions.hpp>
#include <alps/alea/mcdata.hpp>
#include <alps/alea/mcanalyze.hpp>
#include <boost/program_options.hpp>
#include <iostream>

template <class TimeseriesType>
struct return_type {
typedef typename alps::average_type< typename TimeseriesType::value_type >::type type;
};

namespace po = boost::program_options;
typedef std::vector<std::string>::const_iterator const_string_iterator;

template <class X = int>
class impl_calculation {
public:

  template <class TimeseriesType>
  typename return_type<TimeseriesType>::type calculate (TimeseriesType&);

  impl_calculation(std::string name, std::string save_path): _name(name), _save_path(save_path) {}

  int execute (int ac, char* av[]) {
    
    using alps::numeric::operator<<;

    try {
      std::vector<std::string> file_names;
      std::vector<std::string> variable_names;
    
      po::options_description visible("Available options");
      visible.add_options()
        ("help,h", "produce help message")
        ("verbose,v", "print detailed information")
        ("write,w", "write the result(s) back into the file(s)")
        ("name,n", po::value< std::vector<std::string> >(&variable_names), "variable name, can be specified multiple times [default: all variables]")
        ("path,p", po::value<std::string>(&path)->default_value("/simulation/results"), "hdf5-path where the data is stored [default: \"/simulation/results\"]")
        ;
    
      po::options_description hidden("hidden options");
      hidden.add_options()
        ("input-files", po::value< std::vector<std::string> >(&file_names), "input file(s)")
        ;
    
      po::options_description cmdline_options;
      cmdline_options.add(visible).add(hidden);
    
      po::positional_options_description p;
      p.add("input-files", -1);
    
      po::variables_map vm;
      store(po::command_line_parser(ac, av).
        options(cmdline_options).positional(p).run(), vm);
      notify(vm);
    
      if (vm.count("help")) {
        std::cout << visible << std::endl;
        return 0;
      }
    
      if (!vm.count("input-files")) {
        std::cerr << visible << std::endl;
        return 1;
      }

    
      for (const_string_iterator file_it = file_names.begin(); file_it != file_names.end(); ++file_it) {
        alps::hdf5::archive ar(*file_it, "a");
        if (!vm.count("name")) {
          variable_names = ar.list_children(path);
          if (vm.count("verbose")) std::cout << "Variables in file " << *file_it << ":  " << variable_names << std::endl;
        }
        for (const_string_iterator variable_it = variable_names.begin(); variable_it != variable_names.end(); ++variable_it) {

          ar.set_context(path + "/" + ar.encode_segment(*variable_it) +"/");
          if (ar.dimensions("timeseries/data") == 1 ) {
            alps::alea::mcdata<double> data;
            ar >> alps::make_pvp("", data);
            post_calc(calculate(data), ar, *file_it, *variable_it, vm);
          } else if (ar.dimensions("timeseries/data") == 2 ) {
            alps::alea::mcdata<std::vector<double> > data;
            ar >> alps::make_pvp("", data);
            post_calc(calculate(data), ar, *file_it, *variable_it, vm);
          } 
        } // ending for variable_it
      } // ending for file_it
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
      return 1;
    }
   std::cout << "Done\n";
    return 0;
  }


private:

  template <class T>
  void post_calc (const T& result, alps::hdf5::archive& ar, const std::string& file_name, const std::string& variable_name, const po::variables_map& vm) {
    using alps::numeric::operator<<;
    if (vm.count("verbose")) {
      std::cout << "The " << _name << " of variable " << variable_name << " in file " << file_name << " is: " << result << std::endl;
    }
    if (vm.count("write")) {
      ar << alps::make_pvp(path + "/" + ar.encode_segment(variable_name) +"/" + _save_path, result);
    }
  }

  std::string _name;
  std::string _save_path;
  std::string path;

};


#endif
