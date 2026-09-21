// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
#include <alps/solvers.hpp>

NB_MODULE(maxent_c, module) {
  module.def("AnalyticContinuation", [](nanobind::dict const& values) {
    pyalps::scoped_signal_handlers signal_handlers;
    auto parameters = pyalps::params_from_dict(values);
    auto output_file = static_cast<std::string>(parameters["BASENAME"] | "results") + ".out.h5";
    alps::solvers::maxent(parameters, output_file);
  });
}
