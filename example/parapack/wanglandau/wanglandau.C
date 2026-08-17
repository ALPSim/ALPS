/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "wanglandau.h"

PARAPACK_SET_VERSION("ALPS/parapack example program: Wang-Landau algorithm");

PARAPACK_REGISTER_WORKER(wanglandau_worker<alps::wanglandau::learn>, "wanglandau learn");
PARAPACK_REGISTER_EVALUATOR(alps::wanglandau_evaluator<alps::wanglandau::learn>, "wanglandau learn");

PARAPACK_REGISTER_WORKER(wanglandau_worker<alps::wanglandau::measure>, "wanglandau measure");
PARAPACK_REGISTER_EVALUATOR(alps::wanglandau_evaluator<alps::wanglandau::measure>, "wanglandau measure");

PARAPACK_REGISTER_WORKER(wanglandau_reweight_worker, "wanglandau reweight");
PARAPACK_REGISTER_EVALUATOR(wanglandau_reweight_evaluator, "wanglandau reweight");
