/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2013 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/parapack.h>
#include <boost/filesystem/operations.hpp>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Error: " << argv[0] << " master_file\n";
    std::exit(-1);
  }

  alps::oxstream os(std::cout);
  os << alps::header("UTF-8");

  boost::filesystem::path file(argv[1]);
  boost::filesystem::path basedir = file.parent_path();
  std::string file_in_str, file_out_str;

  alps::parapack::load_filename(file, file_in_str, file_out_str);

  typedef std::pair<std::string, std::string> version_t;
  std::vector<version_t> versions;
  alps::parapack::load_version(file, versions);

  boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str), basedir);
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str), basedir);
  std::string simname;
  std::vector<alps::task> tasks;

  alps::parapack::load_tasks(file_in, file_out, basedir, simname, tasks);

  os << alps::start_tag("ARCHIVE")
     << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
     << alps::attribute("xsi:noNamespaceSchemaLocation", "http://xml.comp-phys.org/2008/6/archive.xsd");
  if (simname.size())
    os << alps::attribute("name", simname);

  BOOST_FOREACH(version_t const& v, versions) {
    os << alps::start_tag("VERSION")
       << alps::attribute("type", v.first)
       << alps::attribute("string", v.second)
       << alps::end_tag("VERSION");
  }

  BOOST_FOREACH(alps::task& t, tasks) {
    t.load();
    t.write_xml_archive(os);
  }

  os << alps::end_tag("ARCHIVE");
}
