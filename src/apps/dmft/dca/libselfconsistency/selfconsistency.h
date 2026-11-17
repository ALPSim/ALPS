/*****************************************************************************
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#pragma once
#include "alps/params.hpp"
#include "impurity_solver.h"

void define_selfconsistency_parameters(alps::params &p);

/// set up the cluster transformer and solver and run the DCA self-consistency loop until convergence
void perform_cluster_loop(const std::string &param_file_name, const alps::params &parms);
