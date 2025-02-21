/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/clone.h>

int main(int argc, char** argv) {
  alps::oxstream os(std::cout);
  os << alps::start_tag("CLONES");
  for (int i = 1; i < argc; ++i) {
    alps::hdf5::archive ar(argv[i]);
    alps::Parameters params;
    alps::clone_info info;
    ar["/parameters"] >> params;
    ar["/log/alps"] >> info;
    std::vector<alps::ObservableSet> obs;
    alps::load_observable(ar, info.clone_id(), obs);
    os << alps::start_tag("CLONE")
       << alps::attribute("dumpfile", argv[i])
       << params;
    BOOST_FOREACH(alps::ObservableSet const& m, obs) m.write_xml(os);
    info.write_xml(os);
    os << alps::end_tag("CLONE");
  }
  os << alps::end_tag("CLONES");
}
