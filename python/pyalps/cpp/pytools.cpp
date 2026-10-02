// Copyright (C) 1994-2009 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
//               Matthias Troyer <troyer@itp.phys.ethz.ch>,
//               Bela Bauer <bauerb@itp.phys.ethz.ch>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/scheduler/convert.h>
#include <alps/utility/encode.hpp>
#include <alps/random.h>
#include <alps/parser/xslt_path.h>
#include <alps/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
namespace nb = nanobind;
typedef boost::variate_generator<boost::mt19937, boost::uniform_01<double> > random_01;
class WrappedRNG : public random_01
{
public:
    WrappedRNG(int seed = 0)
        : random_01(boost::mt19937(seed), boost::uniform_01<double>())
    {
    }
};
NB_MODULE(pytools_c, m) {
    m.doc() = "ALPS tools bindings (nanobind)";
    m.def("convert2xml",
          &alps::convert2xml,
          "Convert an ALPS file to XML. Returns the path to the XML file.");
    m.def("hdf5_name_encode",
          &alps::hdf5_name_encode,
          "Escape a string for use inside an HDF5 path name.");
    m.def("hdf5_name_decode",
          &alps::hdf5_name_decode,
          "Un-escape a string taken from an HDF5 path name.");
    m.def("search_xml_library_path",
          &alps::search_xml_library_path,
          "Resolve an ALPS library XML / XSL file to its full path.");
    nb::class_<WrappedRNG>(m, "rng",
                           "Mersenne-Twister uniform random number generator in [0, 1).")
        .def(nb::init<int>(), nb::arg("seed") = 0)
        .def("__deepcopy__",
             [](WrappedRNG const & self, nb::handle /*memo*/) {
                 return WrappedRNG(self);
             },
             "Return a fresh copy of the RNG carrying the same state.")
        .def("__call__",
             static_cast<WrappedRNG::result_type (WrappedRNG::*)()>(
                 &WrappedRNG::operator()),
             "Return a uniform random number in [0, 1).");
}
