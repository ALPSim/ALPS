/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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


#include <alps/hdf5/archive.hpp>
#include <alps/ngs/accumulator.hpp>
#include <alps/ngs/boost_python.hpp>

#include <boost/python.hpp>
#include <boost/optional.hpp>

#include <string>
#include <sstream>

namespace alps {
	namespace accumulator {
		namespace python {

			class object_wrapper {
				public:

					object_wrapper(boost::python::object arg)
						: obj(arg)
					{}

					boost::python::object operator+(boost::python::object arg) {
						return boost::python::call_method<boost::python::object>(obj.get().ptr(), "__add__", arg);
					}

				private:
					boost::optional<boost::python::object> obj;
			};

	        typedef impl::Accumulator<python::object_wrapper, count_tag, impl::AccumulatorBase<python::object_wrapper> > count_accumulator;

			void magic_call(count_accumulator & self, boost::python::object arg) { self(object_wrapper(arg)); }

			template<typename T> std::string magic_str(T & self) { 
				std::stringstream ss; 
				self.print(ss); 
				return ss.str(); 
			}

			count_accumulator::result_type result(count_accumulator & self) { return count_accumulator::result_type(self); }

			count_accumulator::result_type add_result(count_accumulator::result_type self, count_accumulator::result_type const & arg) { self += arg; return self; }
			count_accumulator::result_type add_double(count_accumulator::result_type self, double arg) { self += arg; return self; }

			count_accumulator::result_type sub_result(count_accumulator::result_type self, count_accumulator::result_type const & arg) { self -= arg; return self; }
			count_accumulator::result_type sub_double(count_accumulator::result_type self, double arg) { self -= arg; return self; }
			count_accumulator::result_type rsub_double(count_accumulator::result_type self, double arg) { self.negate(); self += arg; return self; }

			count_accumulator::result_type mul_result(count_accumulator::result_type self, count_accumulator::result_type const & arg) { self *= arg; return self; }
			count_accumulator::result_type mul_double(count_accumulator::result_type self, double arg) { self *= arg; return self; }

			count_accumulator::result_type div_result(count_accumulator::result_type self, count_accumulator::result_type const & arg) { self /= arg; return self; }
			count_accumulator::result_type div_double(count_accumulator::result_type self, double arg) { self /= arg; return self; }
			count_accumulator::result_type rdiv_double(count_accumulator::result_type self, double arg) { self.inverse(); self *= arg; return self; }

			count_accumulator sin(count_accumulator::result_type self) { self.sin(); return self; }
			count_accumulator cos(count_accumulator::result_type self) { self.cos(); return self; }
			count_accumulator tan(count_accumulator::result_type self) { self.tan(); return self; }
			count_accumulator sinh(count_accumulator::result_type self) { self.sinh(); return self; }
			count_accumulator cosh(count_accumulator::result_type self) { self.cosh(); return self; }
			count_accumulator tanh(count_accumulator::result_type self) { self.tanh(); return self; }
			count_accumulator asin(count_accumulator::result_type self) { self.asin(); return self; }
			count_accumulator acos(count_accumulator::result_type self) { self.acos(); return self; }
			count_accumulator atan(count_accumulator::result_type self) { self.atan(); return self; }
			count_accumulator abs(count_accumulator::result_type self) { self.abs(); return self; }
			count_accumulator sqrt(count_accumulator::result_type self) { self.sqrt(); return self; }
			count_accumulator log(count_accumulator::result_type self) { self.log(); return self; }
			count_accumulator sq(count_accumulator::result_type self) { self.sq(); return self; }
			count_accumulator cb(count_accumulator::result_type self) { self.cb(); return self; }
			count_accumulator cbrt(count_accumulator::result_type self) { self.cbrt(); return self; }

		}
	}
}

BOOST_PYTHON_MODULE(pyngsaccumulator_c) {

	using namespace boost::python;

	typedef alps::accumulator::python::count_accumulator count_accumulator;
    class_<count_accumulator>(
        "count_accumulator",
		init<>()
    )
        .def("__call__", &alps::accumulator::python::magic_call)
        .def("__str__", &alps::accumulator::python::magic_str<count_accumulator>)

        .def("save", &count_accumulator::save)
        .def("load", &count_accumulator::load)

        .def("reset", &count_accumulator::reset)

        .def("result", &alps::accumulator::python::result)
    ;

    typedef count_accumulator::result_type count_result_type; 
    class_<count_result_type>(
        "count_result",
		init<>()
    )
        .def("__str__", &alps::accumulator::python::magic_str<count_result_type>)

        .def("save", &count_result_type::save)
        .def("load", &count_result_type::load)

        .def("reset", &count_result_type::reset)

        .def(self += count_result_type())
        .def(self += double())
        .def("__add__", &alps::accumulator::python::add_result)
        .def("__add__", &alps::accumulator::python::add_double)
        .def("__radd__", &alps::accumulator::python::add_double)

        .def(self -= count_result_type())
        .def(self -= double())
        .def("__sub__", &alps::accumulator::python::sub_result)
        .def("__sub__", &alps::accumulator::python::sub_double)
        .def("__rsub__", &alps::accumulator::python::rsub_double)

        .def(self *= count_result_type())
        .def(self *= double())
        .def("__mul__", &alps::accumulator::python::mul_result)
        .def("__mul__", &alps::accumulator::python::mul_double)
        .def("__rmul__", &alps::accumulator::python::mul_double)

        .def(self /= count_result_type())
        .def(self /= double())
        .def("__div__", &alps::accumulator::python::div_result)
        .def("__div__", &alps::accumulator::python::div_double)
        .def("__rdiv__", &alps::accumulator::python::rdiv_double)

        .def("sin", &alps::accumulator::python::sin)
        .def("cos", &alps::accumulator::python::cos)
        .def("tan", &alps::accumulator::python::tan)
        .def("sinh", &alps::accumulator::python::sinh)
        .def("cosh", &alps::accumulator::python::cosh)
        .def("tanh", &alps::accumulator::python::tanh)
        .def("asin", &alps::accumulator::python::asin)
        .def("acos", &alps::accumulator::python::acos)
        .def("atan", &alps::accumulator::python::atan)
        .def("abs", &alps::accumulator::python::abs)
        .def("sqrt", &alps::accumulator::python::sqrt)
        .def("log", &alps::accumulator::python::log)
        .def("sq", &alps::accumulator::python::sq)
        .def("cb", &alps::accumulator::python::cb)
        .def("cbrt", &alps::accumulator::python::cbrt)
    ;
}
