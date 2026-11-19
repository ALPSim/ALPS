/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#pragma once

namespace alps { namespace alea {

/** Outcome of the check whether an error estimate has stopped growing */
enum error_convergence { CONVERGED = 0, MAYBE_CONVERGED = 1, NOT_CONVERGED = 2 };

}}
