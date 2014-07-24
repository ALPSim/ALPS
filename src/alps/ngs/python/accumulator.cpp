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


#include <alps/hdf5/python.hpp>
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
					object_wrapper() {}
					object_wrapper(boost::python::object arg): obj(arg) {}
					object_wrapper(boost::uint64_t arg): obj(arg) {}

					boost::python::object & get() { return obj; }
					boost::python::object const &get() const { return obj; }

					void print(std::ostream & os) const {
						os << boost::python::call_method<std::string>(obj.ptr(), "__str__");
					}

					#define ALPS_ACCUMULATOR_PYTHON_OPERATOR(cxxiop, pyipo, cxxop, pyop)							\
						object_wrapper & cxxiop (object_wrapper arg) {												\
							boost::python::call_method<boost::python::object>(obj.ptr(), pyipo , arg.obj);			\
							return *this;																			\
						}																							\
						boost::python::object cxxop (object_wrapper arg) {											\
							return boost::python::call_method<boost::python::object>(obj.ptr(), pyop, arg.obj);		\
						}

					ALPS_ACCUMULATOR_PYTHON_OPERATOR(operator+=, "__iadd__", operator+, "__add__")
					ALPS_ACCUMULATOR_PYTHON_OPERATOR(operator-=, "__isub__", operator-, "__sub__")
					ALPS_ACCUMULATOR_PYTHON_OPERATOR(operator*=, "__imul__", operator*, "__mul__")
					ALPS_ACCUMULATOR_PYTHON_OPERATOR(operator/=, "__idiv__", operator/, "__div__")
					#undef ALPS_ACCUMULATOR_PYTHON_OPERATOR

				private:
					boost::python::object obj;
			};

			inline std::ostream & operator<<(std::ostream & os, object_wrapper const & arg) {
				arg.print(os);
				return os;
			}

			template<typename T> void magic_call(T & self, boost::python::object arg) { self(object_wrapper(arg)); }

			template<typename T> std::string magic_str(T & self) { 
				std::stringstream ss; 
				self.print(ss); 
				return ss.str(); 
			}

			template<typename T> typename T::result_type result(T & self) { return typename T::result_type(self); }

			template<typename T> typename T::result_type add_result(typename T::result_type self, typename T::result_type const & arg) { self += arg; return self; }
			template<typename T> typename T::result_type add_double(typename T::result_type self, double arg) { self += arg; return self; }

			template<typename T> typename T::result_type sub_result(typename T::result_type self, typename T::result_type const & arg) { self -= arg; return self; }
			template<typename T> typename T::result_type sub_double(typename T::result_type self, double arg) { self -= arg; return self; }
			template<typename T> typename T::result_type rsub_double(typename T::result_type self, double arg) { self.negate(); self += arg; return self; }

			template<typename T> typename T::result_type mul_result(typename T::result_type self, typename T::result_type const & arg) { self *= arg; return self; }
			template<typename T> typename T::result_type mul_double(typename T::result_type self, double arg) { self *= arg; return self; }

			template<typename T> typename T::result_type div_result(typename T::result_type self, typename T::result_type const & arg) { self /= arg; return self; }
			template<typename T> typename T::result_type div_double(typename T::result_type self, double arg) { self /= arg; return self; }
			template<typename T> typename T::result_type rdiv_double(typename T::result_type self, double arg) { self.inverse(); self *= arg; return self; }

			template<typename T> T sin(typename T::result_type self) { self.sin(); return self; }
			template<typename T> T cos(typename T::result_type self) { self.cos(); return self; }
			template<typename T> T tan(typename T::result_type self) { self.tan(); return self; }
			template<typename T> T sinh(typename T::result_type self) { self.sinh(); return self; }
			template<typename T> T cosh(typename T::result_type self) { self.cosh(); return self; }
			template<typename T> T tanh(typename T::result_type self) { self.tanh(); return self; }
			template<typename T> T asin(typename T::result_type self) { self.asin(); return self; }
			template<typename T> T acos(typename T::result_type self) { self.acos(); return self; }
			template<typename T> T atan(typename T::result_type self) { self.atan(); return self; }
			template<typename T> T abs(typename T::result_type self) { self.abs(); return self; }
			template<typename T> T sqrt(typename T::result_type self) { self.sqrt(); return self; }
			template<typename T> T log(typename T::result_type self) { self.log(); return self; }
			template<typename T> T sq(typename T::result_type self) { self.sq(); return self; }
			template<typename T> T cb(typename T::result_type self) { self.cb(); return self; }
			template<typename T> T cbrt(typename T::result_type self) { self.cbrt(); return self; }

		}
	}

	namespace hdf5 {

        template<> struct scalar_type<alps::accumulator::python::object_wrapper> {
            typedef alps::accumulator::python::object_wrapper type;
        };

		namespace detail {

            template<> struct is_vectorizable<alps::accumulator::python::object_wrapper> {
                static bool apply(alps::accumulator::python::object_wrapper const & value) {
                	return is_vectorizable<boost::python::object>::apply(value.get());
                }
            };

            template<> struct get_extent<alps::accumulator::python::object_wrapper> {
                static std::vector<std::size_t> apply(alps::accumulator::python::object_wrapper const & value) {
                	return get_extent<boost::python::object>::apply(value.get());
                }
            };

            template<> struct set_extent<alps::accumulator::python::object_wrapper> {
                static void apply(alps::accumulator::python::object_wrapper & value, std::vector<std::size_t> const & extent) {
                	set_extent<boost::python::object>::apply(value.get(), extent);
                }
            };
        }

        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , alps::accumulator::python::object_wrapper const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
	        save(ar, path, value.get(), size, chunk, offset);
        }

        ALPS_DECL void load(
              archive & ar
            , std::string const & path
            , alps::accumulator::python::object_wrapper & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
        	load(ar, path, value.get(), chunk, offset);
        }
	}
}

BOOST_PYTHON_MODULE(pyngsaccumulator_c) {

	using namespace boost::python;
	using namespace alps::accumulator::impl;

	#define ALPS_ACCUMULATOR_COMMON(accumulator_type)									\
        .def("__str__", &alps::accumulator::python::magic_str< accumulator_type >)		\
        .def("save", & accumulator_type ::save)											\
        .def("load", & accumulator_type ::load)											\
        .def("reset", & accumulator_type ::reset)

	#define ALPS_RESULT_COMMON(accumulator_type)										\
		.def(self += count_result_type())												\
        .def(self += double())															\
        .def("__add__", &alps::accumulator::python::add_result< accumulator_type >)		\
        .def("__add__", &alps::accumulator::python::add_double< accumulator_type >)		\
        .def("__radd__", &alps::accumulator::python::add_double< accumulator_type >)	\
        .def(self -= count_result_type())												\
        .def(self -= double())															\
        .def("__sub__", &alps::accumulator::python::sub_result< accumulator_type >)		\
        .def("__sub__", &alps::accumulator::python::sub_double< accumulator_type >)		\
        .def("__rsub__", &alps::accumulator::python::rsub_double< accumulator_type >)	\
        .def(self *= count_result_type())												\
        .def(self *= double())															\
        .def("__mul__", &alps::accumulator::python::mul_result< accumulator_type >)		\
        .def("__mul__", &alps::accumulator::python::mul_double< accumulator_type >)		\
        .def("__rmul__", &alps::accumulator::python::mul_double< accumulator_type >)	\
        .def(self /= count_result_type())												\
        .def(self /= double())															\
        .def("__div__", &alps::accumulator::python::div_result< accumulator_type >)		\
        .def("__div__", &alps::accumulator::python::div_double< accumulator_type >)		\
        .def("__rdiv__", &alps::accumulator::python::rdiv_double< accumulator_type >)	\
        .def("sin", &alps::accumulator::python::sin< accumulator_type >)				\
        .def("cos", &alps::accumulator::python::cos< accumulator_type >)				\
        .def("tan", &alps::accumulator::python::tan< accumulator_type >)				\
        .def("sinh", &alps::accumulator::python::sinh< accumulator_type >)				\
        .def("cosh", &alps::accumulator::python::cosh< accumulator_type >)				\
        .def("tanh", &alps::accumulator::python::tanh< accumulator_type >)				\
        .def("asin", &alps::accumulator::python::asin< accumulator_type >)				\
        .def("acos", &alps::accumulator::python::acos< accumulator_type >)				\
        .def("atan", &alps::accumulator::python::atan< accumulator_type >)				\
        .def("abs", &alps::accumulator::python::abs< accumulator_type >)				\
        .def("sqrt", &alps::accumulator::python::sqrt< accumulator_type >)				\
        .def("log", &alps::accumulator::python::log< accumulator_type >)				\
        .def("sq", &alps::accumulator::python::sq< accumulator_type >)					\
        .def("cb", &alps::accumulator::python::cb< accumulator_type >)					\
        .def("cbrt", &alps::accumulator::python::cbrt< accumulator_type >)

	typedef alps::accumulator::python::object_wrapper python_object;

	typedef Accumulator<python_object, alps::accumulator::count_tag, AccumulatorBase<python_object> > count_accumulator;
    class_<count_accumulator>(
        "count_accumulator",
		init<>()
    )
        .def("__call__", &alps::accumulator::python::magic_call<count_accumulator>)
    	ALPS_ACCUMULATOR_COMMON(count_accumulator)
        .def("result", &alps::accumulator::python::result<count_accumulator>)

        .def("count", &count_accumulator::count)
    ;

    typedef count_accumulator::result_type count_result_type; 
    class_<count_result_type>(
        "count_result",
		init<>()
    )
    	ALPS_ACCUMULATOR_COMMON(count_result_type)

        .def("count", &count_accumulator::count)

        ALPS_RESULT_COMMON(count_accumulator)
    ;

	typedef Accumulator<python_object, alps::accumulator::mean_tag, count_accumulator> mean_accumulator;
    class_<mean_accumulator>(
        "mean_accumulator",
		init<>()
    )
        .def("__call__", &alps::accumulator::python::magic_call<mean_accumulator>)
    	ALPS_ACCUMULATOR_COMMON(mean_accumulator)
        .def("result", &alps::accumulator::python::result<mean_accumulator>)

        .def("count", &mean_accumulator::count)
        .def("mean", &mean_accumulator::count) // TODO: make function
    ;

    typedef mean_accumulator::result_type mean_result_type; 
    class_<mean_result_type>(
        "mean_result",
		init<>()
    )
    	ALPS_ACCUMULATOR_COMMON(mean_result_type)

        .def("count", &mean_accumulator::count)
        .def("mean", &mean_accumulator::count) // TODO: make function

        // ALPS_RESULT_COMMON(mean_accumulator) // TODO: implement
    ;    
}
