/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define PY_ARRAY_UNIQUE_SYMBOL pyising_PyArrayHandle

#include "ising.hpp"

#include <alps/ngs/boost_python.hpp>

#include <alps/python/make_copy.hpp>

#include <boost/bind.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/return_internal_reference.hpp>

BOOST_PYTHON_MODULE(ngsising_c) {
	boost::python::class_<ising_sim, boost::noncopyable, boost::python::bases<alps::mcbase_ng> >(
	      "sim" ,
	      boost::python::init<ising_sim::parameters_type const &, boost::python::optional<std::size_t> >()
	)
	    .add_property("params", boost::python::make_function(
	        static_cast<ising_sim::parameters_type &(ising_sim::*)()>(&ising_sim::get_params), boost::python::return_internal_reference<>()
	     ))
	    .def("run", static_cast<bool(ising_sim::*)(boost::python::object)>(&ising_sim::run))
	    .def("random", &ising_sim::get_random)
	    .def("save", static_cast<void(ising_sim::*)(alps::hdf5::archive &) const>(&ising_sim::save))
	    .def("load", static_cast<void(ising_sim::*)(alps::hdf5::archive &)>(&ising_sim::load))
	;
}