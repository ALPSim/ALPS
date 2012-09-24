/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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
