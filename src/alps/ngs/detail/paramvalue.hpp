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

#ifndef ALPS_NGS_DETAIL_PARAMVALUE_HPP
#define ALPS_NGS_DETAIL_PARAMVALUE_HPP

#include <alps/ngs/hdf5.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/detail/paramvalue_reader.hpp>

#if defined(ALPS_HAVE_PYTHON)
	#include <alps/ngs/boost_python.hpp>
#endif

#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp> 
#include <boost/serialization/complex.hpp> 
#include <boost/serialization/split_member.hpp>

#include <string>
#include <complex>
#include <ostream>
#include <stdexcept>

#define ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(CALLBACK)					\
	CALLBACK(double)																\
	CALLBACK(int)																	\
	CALLBACK(bool)																	\
	CALLBACK(std::string)															\
	CALLBACK(std::complex<double>)													\
	CALLBACK(std::vector<double>)													\
	CALLBACK(std::vector<int>)														\
	CALLBACK(std::vector<std::string>)												\
	CALLBACK(std::vector<std::complex<double> >)

#if defined(ALPS_HAVE_PYTHON)
	#define ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(CALLBACK)							\
		ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(CALLBACK)					\
		CALLBACK(boost::python::object)
#else
	#define ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(CALLBACK)							\
		ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(CALLBACK)
#endif

namespace alps {

	namespace detail {

		class paramvalue;

	}

	template<typename T> class extract;

	namespace detail {

		template<class Archive> struct paramvalue_serializer 
			: public boost::static_visitor<> 
		{
			public:

				paramvalue_serializer(Archive & a)
					: ar(a)
				{}

				template <typename U> void operator()(U & v) const {
					std::string type(typeid(v).name());
					ar
						<< type
						<< v;
				}

			private:

				Archive & ar;
        };

		#define ALPS_NGS_PARAMVALUE_VARIANT_TYPE(T)	T,
		typedef boost::mpl::pop_back<boost::mpl::vector<
			ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(
				ALPS_NGS_PARAMVALUE_VARIANT_TYPE
			) void
		>::type >::type paramvalue_types;
		#undef ALPS_NGS_PARAMVALUE_VARIANT_TYPE
        typedef boost::make_variant_over<paramvalue_types>::type paramvalue_base;

		class paramvalue : public paramvalue_base {

			public:

				paramvalue() {}

				paramvalue(paramvalue const & v)
					: paramvalue_base(static_cast<paramvalue_base const &>(v))
				{}

				template<typename T> T cast() const {
					return extract<T>(*this);
				}

				#define ALPS_NGS_PARAMVALUE_MEMBER_DECL(T)							\
					paramvalue( T const & v) : paramvalue_base(v) {}				\
					operator T () const;											\
					paramvalue & operator=( T const &);
				ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_MEMBER_DECL)
				#undef ALPS_NGS_PARAMVALUE_MEMBER_DECL

				void save(hdf5::archive &) const;
				void load(hdf5::archive &);
				
			private:
			
				friend class boost::serialization::access;

				template<class Archive> void save(
					Archive & ar, const unsigned int version
				) const {
					paramvalue_serializer<Archive> visitor(ar);
					boost::apply_visitor(visitor, *this);
				}

				template<class Archive> void load(
					Archive & ar, const unsigned int
				) {
					std::string type;
					ar >> type;
					if (false);
					#define ALPS_NGS_PARAMVALUE_LOAD(T)								\
						else if (type == typeid( T ).name()) {						\
							T value;												\
							ar >> value;											\
							operator=(value);										\
						}
					ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_LOAD)
					#undef ALPS_NGS_PARAMVALUE_LOAD
					else
						throw std::runtime_error("unknown type" + ALPS_STACKTRACE);
				}

				BOOST_SERIALIZATION_SPLIT_MEMBER()
		};

		std::ostream & operator<<(std::ostream & os, paramvalue const & arg);

		template<typename T> T extract_impl (paramvalue const & arg, T) {
			paramvalue_reader< T > visitor;
			boost::apply_visitor(visitor, arg);
			return visitor.get_value();
		}

		template<typename T> struct convert_hook<T, paramvalue> {
			static inline std::complex<T> apply(paramvalue const & arg) {
				return extract<T>(arg);
			}
		};

    }

	template<typename T> class extract {

		public:

			extract(detail::paramvalue const & arg) {
				boost::apply_visitor(visitor, arg);
			}

			operator T const & () {
				return visitor.get_value();
			}

		private:

			detail::paramvalue_reader< T > visitor;
	};
}

#endif
