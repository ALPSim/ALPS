/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parser/xmlstream.h>
#include <alps/parser/xslt_path.h>
#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <string>

namespace po = boost::program_options;
using namespace alps;

int main(int argc, char** argv) {

  po::options_description desc("Required/allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("name,n", po::value<std::string>(), "name of simulation")
    ("version,v", po::value<std::string>(), "version information of application")
    ("xaxis,x", po::value<std::string>(), "name of xaxis, or parameter [required]")
    ("yaxis,y", po::value<std::string>(), "name of yaxis, or observable [required]")
    ("with-error,e", "specify if observable has error bar")
    ("input-file,i", po::value<std::string>(), "input text file [required]");
  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cerr << desc << std::endl;
    std::exit(0);
  }

  if (!vm.count("input-file") || !vm.count("xaxis") || !vm.count("yaxis")) {
    std::cerr << "txt2archive: invalid options\n" << desc << std::endl;
    std::exit(-1);
  }

  // read from file
  std::ifstream is(vm["input-file"].as<std::string>().c_str());
  bool with_error = vm.count("with-error");

  // output to Archive XML
  std::string xaxis = vm["xaxis"].as<std::string>();
  std::string yaxis = vm["xaxis"].as<std::string>();

  oxstream os(std::cout);
  os << header("UTF-8");

  os << start_tag("ARCHIVE")
     << xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
     << attribute("xsi:noNamespaceSchemaLocation",
                        "http://xml.comp-phys.org/2008/6/archive.xsd");
  if (vm.count("name"))
    os << attribute("name", vm["name"].as<std::string>());

  if (vm.count("version"))
    os << start_tag("VERSION")
       << attribute("type", "application")
       << attribute("string", vm["version"].as<std::string>())
       << end_tag("VERSION");

  std::string str;
  int id = 0;
  while (std::getline(is, str) && str.size()) {
    std::istringstream iss(str);
    double x, y, dy;
    if (with_error) {
      iss >> x >> y >> dy;
    } else {
      iss >> x >> y;
    }
    os << start_tag("SIMULATION");
    os << attribute("id", ++id);

    os << start_tag("PARAMETERS")
       << start_tag("PARAMETER") << no_linebreak << attribute("name", xaxis) << x
       << end_tag("PARAMETER")
       << end_tag("PARAMETERS");

    os << start_tag("AVERAGES")
       << start_tag("SCALAR_AVERAGE")
       << attribute("name", yaxis)
       << start_tag("MEAN") << no_linebreak << y << end_tag("MEAN");
    if (with_error)
      os << start_tag("ERROR") << no_linebreak << dy << end_tag("ERROR");
    os << end_tag("SCALAR_AVERAGE")
       << end_tag("AVERAGES");

    os << end_tag("SIMULATION");
  }


  os << end_tag("ARCHIVE");
}
