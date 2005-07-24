/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#include <alps/xml.h>
#include <alps/parameterlist.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

void convert_params(const std::string& inname, const std::string& basename)
{
  alps::ParameterList list;
  {
    std::ifstream in(inname.c_str());
    in >> list;
  }

  if (list.size() == 0) return;

  std::cout << "Converting parameter file " << inname << " to "
            <<  basename+".in.xml" << std::endl;

  int bits = 31;
  for (int n = 1; n < list.size(); n<<=1, --bits);

  uint32_t baseseed;
  if (list[0].defined("BASESEED")) {
    baseseed = boost::lexical_cast<uint32_t>(list[0]["BASESEED"]);
  } else {
    baseseed =
      boost::posix_time::microsec_clock::local_time().time_of_day().
      total_microseconds() << 24;
  }
  baseseed &= ((1<<30)|((1<<30)-1));

  alps::oxstream out(boost::filesystem::path((basename+".in.xml").c_str(),boost::filesystem::native));
  out << alps::header("UTF-8")
      << alps::stylesheet(alps::xslt_path("job.xsl"))
      << alps::start_tag("JOB")
      << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
      << alps::attribute("xsi:noNamespaceSchemaLocation",
                         "http://xml.comp-phys.org/2003/8/job.xsd")
      << alps::start_tag("OUTPUT")
      << alps::attribute("file", basename+".out.xml")
      << alps::end_tag("OUTPUT");

  for (int i = 0; i < list.size(); ++i) {
    std::string taskname =
      basename+".task"+boost::lexical_cast<std::string,int>(i+1);

    if (!list[i].defined("SEED")) {
      uint32_t seed = baseseed ^ (i << bits);
      list[i]["SEED"] = seed;
    }

    out << alps::start_tag("TASK") << alps::attribute("status","new")
        << alps::start_tag("INPUT")
        << alps::attribute("file", taskname + ".in.xml")
        << alps::end_tag("INPUT")
        << alps::start_tag("OUTPUT")
        << alps::attribute("file", taskname + ".out.xml")
        << alps::end_tag("OUTPUT")
        << alps::end_tag("TASK");

    alps::oxstream task(boost::filesystem::path((taskname+".in.xml").c_str(),boost::filesystem::native));
    task << alps::header("UTF-8")
         << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
    task << alps::start_tag("SIMULATION")
         << alps::xml_namespace("xsi",
                                "http://www.w3.org/2001/XMLSchema-instance")
         << alps::attribute("xsi:noNamespaceSchemaLocation",
                            "http://xml.comp-phys.org/2002/10/QMCXML.xsd");
    task << list[i];
    task << alps::end_tag("SIMULATION");
  }

  out << alps::end_tag("JOB");
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  if (argc<2 || argc>3) {
    std::cerr << "Usage: " << argv[0] << " inputfile [outputbasename]\n";
    std::exit(-1);
  }

  std::string inname = argv[1];
  std::string outbase = argv[argc-1];
  if (inname.size() >= 2 && inname.substr(0, 2) == "./") inname.erase(0, 2);
  if (outbase.size() >= 2 && outbase.substr(0, 2) == "./") outbase.erase(0, 2);

  convert_params(inname, outbase);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif

}
