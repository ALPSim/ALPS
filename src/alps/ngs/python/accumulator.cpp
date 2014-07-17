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


#include <alps/ngs/accumulator.hpp>
#include <alps/ngs/boost_python.hpp>

#include <boost/optional.hpp>
#include <boost/python/object.hpp>

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

			void magic_call(count_accumulator & self, boost::python::object arg) {
				self(object_wrapper(arg));
			}

			std::string magic_str(count_accumulator & self) {
				std::stringstream ss;
				self.print(ss);
				return ss.str();
			}

		}
	}
}

BOOST_PYTHON_MODULE(pyngsaccumulator_c) {

    boost::python::class_<alps::accumulator::python::count_accumulator>(
        "accumulator",
		boost::python::init<>()
    )
        .def("__call__", &alps::accumulator::python::magic_call)
        .def("__str__", &alps::accumulator::python::magic_str)
    ;

}
