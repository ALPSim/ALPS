// Copyright (C) 2010 by Matthias Troyer
//               2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_PYTHON_MAKE_COPY_HPP
#define ALPS_PYTHON_MAKE_COPY_HPP

#include <nanobind/nanobind.h>

namespace alps {
namespace python {

template <class T>
T make_copy(T const & value, nanobind::handle /*memo*/) {
    return value;
}

}  // namespace python
}  // namespace alps

#endif
