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

#include <alps/hdf5/python.hpp>
#include <alps/hdf5/archive.hpp>

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
					template<typename T> object_wrapper(T arg, typename boost::enable_if<boost::is_scalar<T> >::type* = NULL): obj(arg) {}
					object_wrapper(boost::python::object arg): obj(arg) {}

					operator boost::python::object() { return obj; }
					operator const boost::python::object() const { return obj; }

					boost::python::object & get() { return obj; }
					boost::python::object const &get() const { return obj; }

					void print(std::ostream & os) const {
						os << boost::python::call_method<std::string>(obj.ptr(), "__str__");
					}

					#define ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR(cxxiop, cxxop, iop, op)					\
						object_wrapper & cxxiop (object_wrapper const arg) {										\
							if (obj == boost::python::object())														\
								obj = arg.obj;																		\
							else																					\
								obj iop arg.obj; 																	\
							return *this;																			\
						}																							\
						object_wrapper & cxxiop (double arg) {														\
							if (obj == boost::python::object())														\
								obj = boost::python::object(arg);													\
							else																					\
								obj iop boost::python::object(arg); 												\
							return *this;																			\
						}																							\
						object_wrapper cxxop (object_wrapper const arg) const {										\
							return obj op arg.obj;																	\
						}																							\
						object_wrapper cxxop (double arg) const {													\
							return obj op boost::python::object(arg);												\
						}
					ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR(operator+=, operator+, +=, +)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR(operator-=, operator-, -=, -)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR(operator*=, operator*, *=, *)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR(operator/=, operator/, /=, /)
					#undef ALPS_ACCUMULATOR_PYTHON_MEMBER_NUMERIC_OPERATOR

					#define ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(cxxop, op)							\
						bool cxxop (object_wrapper const arg) const {												\
							return obj op arg.obj;																	\
						}
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator==, ==)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator!=, !=)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator<, <)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator<=, <=)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator>, >)
					ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR(operator>=, >=)
					#undef ALPS_ACCUMULATOR_PYTHON_MEMBER_COMARISON_OPERATOR

					object_wrapper operator- () {
						return boost::python::call_method<boost::python::object>(obj.ptr(), "__neg__");
					}

				private:
					boost::python::object obj;
			};

			inline std::ostream & operator<<(std::ostream & os, object_wrapper const & arg) {
				arg.print(os);
				return os;
			}

			#define ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR(cxxop, op)												\
				inline object_wrapper cxxop (double arg1, object_wrapper const & arg2) {							\
					return boost::python::object(arg1) op arg2.get();												\
				}
			ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR(operator+, +)
			ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR(operator-, -)
			ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR(operator*, *)
			ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR(operator/, /)
			#undef ALPS_ACCUMULATOR_FREE_MEMBER_OPERATOR

			template<typename T> void magic_call(T & self, boost::python::object arg) { self(object_wrapper(arg)); }

			template<typename T> std::string magic_str(T & self) { 
				std::stringstream ss; 
				self.print(ss); 
				return ss.str(); 
			}

			#define ALPS_ACCUMULATOR_PYTHON_FUNCTION(name)															\
				object_wrapper name (object_wrapper const & arg) {													\
					boost::python::object np = boost::python::import("numpy");										\
					return boost::python::call_method<boost::python::object>(np.ptr(), #name, arg.get());			\
				}
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(sin)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(cos)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(tan)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(sinh)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(cosh)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(tanh)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(asin)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(acos)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(atan)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(abs)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(sqrt)
			ALPS_ACCUMULATOR_PYTHON_FUNCTION(log)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(sq)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(cb)
			// ALPS_ACCUMULATOR_PYTHON_FUNCTION(cbrt)
			#undef ALPS_ACCUMULATOR_PYTHON_FUNCTION

			template<typename T> typename T::result_type result(T & self) { return typename T::result_type(self); }

			template<typename T> typename T::result_type neg_result(typename T::result_type self) { self.negate(); return self; }

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

			template<typename T> typename T::result_type sin(typename T::result_type self) { self.sin(); return self; }
			template<typename T> typename T::result_type cos(typename T::result_type self) { self.cos(); return self; }
			template<typename T> typename T::result_type tan(typename T::result_type self) { self.tan(); return self; }
			template<typename T> typename T::result_type sinh(typename T::result_type self) { self.sinh(); return self; }
			template<typename T> typename T::result_type cosh(typename T::result_type self) { self.cosh(); return self; }
			template<typename T> typename T::result_type tanh(typename T::result_type self) { self.tanh(); return self; }
			template<typename T> typename T::result_type asin(typename T::result_type self) { self.asin(); return self; }
			template<typename T> typename T::result_type acos(typename T::result_type self) { self.acos(); return self; }
			template<typename T> typename T::result_type atan(typename T::result_type self) { self.atan(); return self; }
			template<typename T> typename T::result_type abs(typename T::result_type self) { self.abs(); return self; }
			template<typename T> typename T::result_type sqrt(typename T::result_type self) { self.sqrt(); return self; }
			template<typename T> typename T::result_type log(typename T::result_type self) { self.log(); return self; }
			// template<typename T> typename T::result_type sq(typename T::result_type self) { self.sq(); return self; }
			// template<typename T> typename T::result_type cb(typename T::result_type self) { self.cb(); return self; }
			// template<typename T> typename T::result_type cbrt(typename T::result_type self) { self.cbrt(); return self; }

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

    namespace ngs {
        namespace numeric {
            template<> struct inf<alps::accumulator::python::object_wrapper> {
                operator alps::accumulator::python::object_wrapper const() {
                	return alps::accumulator::python::object_wrapper(std::numeric_limits<double>::infinity());
                }
            };
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
		.def(self += accumulator_type ::result_type())									\
        .def(self += double())															\
        .def("__neg__", &alps::accumulator::python::neg_result< accumulator_type >)		\
        .def("__add__", &alps::accumulator::python::add_result< accumulator_type >)		\
        .def("__add__", &alps::accumulator::python::add_double< accumulator_type >)		\
        .def("__radd__", &alps::accumulator::python::add_double< accumulator_type >)	\
        .def(self -= accumulator_type ::result_type())									\
        .def(self -= double())															\
        .def("__sub__", &alps::accumulator::python::sub_result< accumulator_type >)		\
        .def("__sub__", &alps::accumulator::python::sub_double< accumulator_type >)		\
        .def("__rsub__", &alps::accumulator::python::rsub_double< accumulator_type >)	\
        .def(self *= accumulator_type ::result_type())									\
        .def(self *= double())															\
        .def("__mul__", &alps::accumulator::python::mul_result< accumulator_type >)		\
        .def("__mul__", &alps::accumulator::python::mul_double< accumulator_type >)		\
        .def("__rmul__", &alps::accumulator::python::mul_double< accumulator_type >)	\
        .def(self /= accumulator_type ::result_type())									\
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
        /*.def("asin", &alps::accumulator::python::asin< accumulator_type >)			\
        .def("acos", &alps::accumulator::python::acos< accumulator_type >)				\
        .def("atan", &alps::accumulator::python::atan< accumulator_type >)*/			\
        .def("abs", &alps::accumulator::python::abs< accumulator_type >)				\
        .def("sqrt", &alps::accumulator::python::sqrt< accumulator_type >)				\
        .def("log", &alps::accumulator::python::log< accumulator_type >)				\
        /*.def("sq", &alps::accumulator::python::sq< accumulator_type >)				\
        .def("cb", &alps::accumulator::python::cb< accumulator_type >)					\
        .def("cbrt", &alps::accumulator::python::cbrt< accumulator_type >)*/

	typedef alps::accumulator::python::object_wrapper python_object;

	typedef Accumulator<python_object, alps::accumulator::count_tag, AccumulatorBase<python_object> > count_accumulator_type;
    class_<count_accumulator_type>("count_accumulator", init<>())
        .def("__call__", &alps::accumulator::python::magic_call<count_accumulator_type>)
    	ALPS_ACCUMULATOR_COMMON(count_accumulator_type)
        .def("result", &alps::accumulator::python::result<count_accumulator_type>)

        .def("count", &count_accumulator_type::count)
    ;

    typedef count_accumulator_type::result_type count_result_type; 
    class_<count_result_type>("count_result", init<>())
    	ALPS_ACCUMULATOR_COMMON(count_result_type)

        .def("count", &count_accumulator_type::count)

        ALPS_RESULT_COMMON(count_accumulator_type)
    ;

	typedef Accumulator<python_object, alps::accumulator::mean_tag, count_accumulator_type> mean_accumulator_type;
    class_<mean_accumulator_type>("mean_accumulator", init<>())
        .def("__call__", &alps::accumulator::python::magic_call<mean_accumulator_type>)
    	ALPS_ACCUMULATOR_COMMON(mean_accumulator_type)
        .def("result", &alps::accumulator::python::result<mean_accumulator_type>)

        .def("count", &mean_accumulator_type::count)
        .def("mean", &mean_accumulator_type::mean)
    ;

    typedef mean_accumulator_type::result_type mean_result_type; 
    class_<mean_result_type>("mean_result", init<>())
    	ALPS_ACCUMULATOR_COMMON(mean_result_type)

        .def("count", &mean_accumulator_type::count)
        .def("mean", &mean_accumulator_type::mean)

        ALPS_RESULT_COMMON(mean_accumulator_type)
    ;

	typedef Accumulator<python_object, alps::accumulator::error_tag, mean_accumulator_type> error_accumulator_type;
    class_<error_accumulator_type>("error_accumulator", init<>())
        .def("__call__", &alps::accumulator::python::magic_call<error_accumulator_type>)
        ALPS_ACCUMULATOR_COMMON(error_accumulator_type)
        .def("result", &alps::accumulator::python::result<error_accumulator_type>)

        .def("count", &error_accumulator_type::count)
        .def("mean", &error_accumulator_type::mean)
        .def("error", &error_accumulator_type::error)
    ;

    typedef error_accumulator_type::result_type error_result_type; 
    class_<error_result_type>("error_result", init<>())
    	ALPS_ACCUMULATOR_COMMON(error_result_type)

        .def("count", &error_accumulator_type::count)
        .def("mean", &error_accumulator_type::mean)
        .def("error", &error_accumulator_type::error)

        ALPS_RESULT_COMMON(error_accumulator_type)
    ;

	typedef Accumulator<python_object, alps::accumulator::binning_analysis_tag, error_accumulator_type> binning_analysis_accumulator_type;
    class_<binning_analysis_accumulator_type>("binning_analysis_accumulator", init<>())
        .def("__call__", &alps::accumulator::python::magic_call<binning_analysis_accumulator_type>)
        ALPS_ACCUMULATOR_COMMON(binning_analysis_accumulator_type)
        .def("result", &alps::accumulator::python::result<binning_analysis_accumulator_type>)

        .def("count", &binning_analysis_accumulator_type::count)
        .def("mean", &binning_analysis_accumulator_type::mean)
        .def("error", &binning_analysis_accumulator_type::error)
    ;

    typedef binning_analysis_accumulator_type::result_type binning_analysis_result_type; 
    class_<binning_analysis_result_type>("binning_analysis_result", init<>())
    	ALPS_ACCUMULATOR_COMMON(binning_analysis_result_type)

        .def("count", &binning_analysis_accumulator_type::count)
        .def("mean", &binning_analysis_accumulator_type::mean)
        .def("error", &binning_analysis_accumulator_type::error)

        ALPS_RESULT_COMMON(binning_analysis_accumulator_type)
    ;    


            // template<typename T> struct observable_type
            //     : public impl::Accumulator<T, max_num_binning_tag, impl::Accumulator<T, binning_analysis_tag, simple_observable_type<T> > >
            // {
            //     typedef typename impl::Accumulator<T, max_num_binning_tag, impl::Accumulator<T, binning_analysis_tag, typename simple_observable_type<T>::accumuator_type> > accumuator_type;
            //     typedef typename impl::Result<T, max_num_binning_tag, impl::Result<T, binning_analysis_tag, typename simple_observable_type<T>::result_type> > result_type;
            //     observable_type(): base_type() {}
            //     template<typename A> observable_type(A const & arg): base_type(arg) {}
            //     private:
            //         typedef impl::Accumulator<T, max_num_binning_tag, impl::Accumulator<T, binning_analysis_tag, simple_observable_type<T> > > base_type;
            // };


}
