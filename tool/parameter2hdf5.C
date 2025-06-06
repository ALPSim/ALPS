/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2015 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*                            Lukas Gamper <gamperl@gmail.com>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include <alps/xml.h>
#include <alps/parser/xslt_path.h>
#include <alps/parameter.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

void convert_params(const std::string& inname, const std::string& basename,
                    const alps::Parameters& params_add)
{
  alps::ParameterList list;
  {
    std::ifstream in(inname.c_str());
    in >> list;
  }
  if (list.size() == 0) {
    alps::Parameters params;
    {
      std::ifstream in(inname.c_str());
      in >> params;
    }
    if (params.size() == 0) {
      std::cerr << "Warning: no parameter set is found.  No HDF5 files will be generated.\n";
      return;
    } else {
      std::cerr << "Info: no parameter set is found.  One task HDF5 will be generated by using default parameters.\n";
      list.push_back(params);
    }
  }

  BOOST_FOREACH(alps::Parameters& params, list) {
    BOOST_FOREACH(alps::Parameter const& padd, params_add) {
      std::string key = padd.key();
      if (params.defined(key)) {
        if (params[key] != padd.value()) {
          std::cerr << "Error: parameter " << key << " is already defined\n";
          boost::throw_exception(std::invalid_argument("convert_params"));
        }
      } else {
        params[key] = padd.value();
      }
    }
  }

  std::cout << "Converting parameter file " << inname << " to "
            <<  basename+".in.xml" << std::endl;

  int bits = 31;
  for (int n = 1; n < list.size(); n<<=1, --bits) {}

  alps::uint32_t baseseed;
  if (list[0].defined("BASESEED")) {
    baseseed = boost::lexical_cast<alps::uint32_t>(list[0]["BASESEED"]);
  } else {
    baseseed =
      boost::posix_time::microsec_clock::local_time().time_of_day().total_microseconds();
    baseseed = ((baseseed << 10) | (baseseed >> 22));
  }

  std::string sim_name = list[0].value_or_default("SIMULATION_NAME", "");
  BOOST_FOREACH(alps::Parameters& p, list) {
    if (p.value_or_default("SIMULATION_NAME", "") != sim_name)
      boost::throw_exception(std::invalid_argument("inconsistent SIMULATION_NAME parameter"));
    if (p.defined("SIMULATION_NAME")) p.erase("SIMULATION_NAME");
  }

  boost::filesystem::path infile((basename+".in.xml").c_str());
  alps::oxstream out(infile);

  // make sure ths stylesheet is there
  alps::copy_stylesheet(infile.remove_filename());



  out << alps::header("UTF-8")
      << alps::stylesheet(alps::xslt_path("ALPS.xsl"))
      << alps::start_tag("JOB")
      << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
      << alps::attribute("xsi:noNamespaceSchemaLocation",
                         "http://xml.comp-phys.org/2003/8/job.xsd");
  if (sim_name != "") out << alps::attribute("name", sim_name);
  out << alps::start_tag("OUTPUT")
      << alps::attribute("file", basename+".out.h5")
      << alps::end_tag("OUTPUT");

  for (int i = 0; i < list.size(); ++i) {
    std::string taskname =
      basename+".task"+boost::lexical_cast<std::string,int>(i+1);

    if (!list[i].defined("SEED")) {
      alps::uint32_t seed = baseseed;
      for (int j = 0; j <= (32/bits); ++j) seed ^= (i << (j * bits));
      seed &= ((1<<30) | ((1<<30)-1));
      list[i]["SEED"] = seed;
    }

    out << alps::start_tag("TASK") << alps::attribute("status","new")
        << alps::start_tag("INPUT")
        << alps::attribute("file", taskname + ".in.h5")
        << alps::end_tag("INPUT")
        << alps::start_tag("OUTPUT")
        << alps::attribute("file", taskname + ".out.h5")
        << alps::end_tag("OUTPUT")
        << alps::end_tag("TASK");

    {
      alps::hdf5::archive ar(boost::filesystem::path(taskname + ".in.h5"), "w");
      ar["/parameters"] << list[i];
    }
  }

  out << alps::end_tag("JOB");
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [-f] inputfile [outputbasename] [[NAME=value]...]\n";
    std::exit(-1);
  }

  bool force = false;
  bool find_in = false;
  bool find_out = false;
  std::string inname;
  std::string outbase;
  alps::Parameters params_add;

  for (int p = 1; p < argc; ++p) {
    std::string arg = argv[p];
    int eqpos = arg.find('=');
    if (eqpos < arg.size()) {
      if ((eqpos > 0) && (eqpos+1 < arg.size())) {
        params_add[arg.substr(0, eqpos)] = arg.substr(eqpos+1);
      } else {
        std::cerr << "Error: invalid parameter \"" << arg << "\"\n";
        std::exit(-1);
      }
    } else if (arg == "-f") {
      force = true;
    } else if (!find_in) {
      inname = arg;
      outbase = arg;
      find_in = true;
    } else if (!find_out) {
      outbase = arg;
    } else {
      std::cerr << "Usage: " << argv[0] << " [-f] inputfile [outputbasename] [[NAME=value]...]\n";
      std::exit(-1);
    }
  }

  if (inname.size() >= 2 && inname.substr(0, 2) == "./") inname.erase(0, 2);
  if (outbase.size() >= 2 && outbase.substr(0, 2) == "./") outbase.erase(0, 2);
  if (outbase.size() >= 2 && outbase.substr(outbase.size()-1, outbase.size()) == "/")
    outbase.erase(outbase.size()-1, outbase.size());
  if (outbase.size() == 0 || outbase == ".")
    outbase = boost::filesystem::path(inname).filename().string();
  if (boost::filesystem::is_directory(outbase))
    outbase += '/' + boost::filesystem::path(inname).filename().string();
  if (!force && boost::filesystem::exists(outbase + ".out.xml")) {
    std::cerr << "Output files (" + outbase + ".out.xml, etc) exist.  "
              << "Please use '-f' option to force replacing input HDF5 files.\n"
              << "Do you really want to overwrite input files? [y/N] ";
    char res;
    std::cin >> res;
    if (res != 'y' && res != 'Y') {
      std::cerr << "Aborted.\n";
      return 255;
    }
  }

  convert_params(inname, outbase, params_add);

  return 0;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif

}
