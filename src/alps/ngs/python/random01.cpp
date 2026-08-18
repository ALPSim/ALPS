/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define PY_ARRAY_UNIQUE_SYMBOL pyngsrandom_PyArrayHandle

#include <alps/ngs/random01.hpp>
#include <alps/ngs/boost_python.hpp>

#include <alps/python/make_copy.hpp>

BOOST_PYTHON_MODULE(pyngsrandom01_c) {

    boost::python::class_<alps::random01>(
        "random01",
        boost::python::init<boost::python::optional<int> >()
    )
        .def("__deepcopy__",  &alps::python::make_copy<alps::random01>)
        .def("__call__", static_cast<alps::random01::result_type(alps::random01::*)()>(&alps::random01::operator()))
        .def("save", &alps::random01::save)
        .def("load", &alps::random01::load)
    ;
}
