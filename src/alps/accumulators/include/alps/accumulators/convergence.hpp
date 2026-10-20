/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

/* $Id: obsvalue.h 3435 2009-11-28 14:45:38Z troyer $ */

#pragma once

#include <alps/config.hpp>
#include <string>
// ALPS: for boost::throw_exception and std::logic_error, so that this header
// compiles on its own.
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {

enum error_convergence {CONVERGED, MAYBE_CONVERGED, NOT_CONVERGED};

template<typename T> inline std::string convergence_to_text(T) {
	boost::throw_exception(std::logic_error("Not Implemented"));
	return std::string();
}

inline std::string convergence_to_text(int c)
{
  return (c==CONVERGED ? "yes" : c==MAYBE_CONVERGED ? "maybe" : c==NOT_CONVERGED ? "no" : "");
}

} // end namespace alps

