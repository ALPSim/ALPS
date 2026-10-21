// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
#include <alps/solvers.hpp>

NB_MODULE(cthyb, module) {
  module.def("solve", [](nanobind::dict const& values) {
    pyalps::scoped_signal_handlers signal_handlers;
    auto parameters = pyalps::params_from_dict(values);
    parameters.define<std::string>("BASENAME", "results", "base name of the output file");
    auto output_file = parameters["BASENAME"].as<std::string>() + ".out.h5";
    alps::solvers::cthyb(parameters, output_file);
  });
}
