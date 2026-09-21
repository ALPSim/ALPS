// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
#include <alps/solvers.hpp>

NB_MODULE(ctint, module) {
  module.def("solve", [](nanobind::dict const& values) {
    pyalps::scoped_signal_handlers signal_handlers;
    auto parameters = pyalps::params_from_dict(values);
    auto output_file = static_cast<std::string>(parameters["BASENAME"] | "results") + ".out.h5";
    alps::solvers::ctint(parameters, output_file);
  });
}
