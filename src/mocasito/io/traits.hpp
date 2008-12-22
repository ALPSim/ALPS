// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <cstddef>
#include <complex>
#include <string>
#include <cstring>
#include <stdexcept>
#ifndef IO_TRAITS
#define IO_TRAITS
namespace mocasito {
	namespace io {
		template<typename T = void> struct type_traits {
			typedef int type;
			static const int value = -1;
		};
		template<typename T = void> struct value_traits {
			typedef typename T::value_type type;
		};
		template<typename T> struct base_type {
			typedef typename base_type<typename T::value_type>::type type;
		};
		#define MOCASITO_IO_TRAITS_CV(CV)													\
			template<typename T> struct value_traits<T CV> {								\
				typedef typename T::value_type CV type;										\
			};																				\
			template<typename T> struct base_type<T CV> {									\
				typedef typename base_type<typename T::value_type>::type CV type;			\
			};
		MOCASITO_IO_TRAITS_CV(const)
		MOCASITO_IO_TRAITS_CV(volatile)
		MOCASITO_IO_TRAITS_CV(const volatile)
		template<typename T> struct has_fixed_size {
			static const bool value = false;
		};
		template <typename T1, typename T2> struct fixed_size_helper {
			static const std::size_t value = 0;
		};
		template<typename T> struct fixed_size
			: public fixed_size_helper<T, T>
		{};
		template<typename T, bool = true> struct get_size_helper {
			static std::size_t apply(T const & v) {
				return fixed_size<T>::value;
			}
		};
		template<typename T> struct get_size_helper<T, false> {
			static std::size_t apply(T const & v) {
				if (has_fixed_size<typename value_traits<T>::type>::value)
					return v.size() * fixed_size<typename value_traits<T>::type>::value;
				else
					throw(std::runtime_error("the type is not supported!"));
			}
		};
		template<typename T, bool = true> struct set_size_helper {
			static void apply(T & v, std::size_t s) {}
		};
		template<typename T> struct set_size_helper<T, false> {
			static void apply(T & v, std::size_t s) {
				if (has_fixed_size<typename value_traits<T>::type>::value)
					return v.resize(s / fixed_size<typename T::value_type>::value);
				else
					throw(std::runtime_error("the type is not supported!"));
			}
		};
		template<typename T> struct get_data {
			static typename value_traits<T>::type * apply(T & v) {
				return &*(v.begin());
			}
			static typename value_traits<T const>::type * apply(T const & v) {
				return &*(v.begin());
			}
			static typename value_traits<T volatile>::type * apply(T volatile & v) {
				return &*(v.begin());
			}
			static typename value_traits<T const volatile>::type * apply(T const volatile & v) {
				return &*(v.begin());
			}
		};
		#define MOCASITO_IO_GET_DATA(T, arg, expr)											\
			template<> struct get_data<T> {													\
				static value_traits<T>::type * apply(T & arg) {								\
					return expr;															\
				}																			\
				static value_traits<T const>::type * apply(T const & arg) {					\
					return expr;															\
				}																			\
			};
		template<typename T> struct get_size {
			static std::size_t apply(T const & v) {
				return get_size_helper<T, has_fixed_size<T>::value>::apply(v);
			}
		};
		template<typename T> struct set_size {
			static void apply(T & v, std::size_t s) {
				set_size_helper<T, has_fixed_size<T>::value>::apply(v, s);
			}
		};
		// Scalar types
		#define MOCASITO_IO_SET_SCALAR_FIXED(T, N)											\
			template<> struct type_traits<T> {												\
				static const int value = N;													\
			};																				\
			template<> struct value_traits<T> {												\
				typedef T type;																\
			};																				\
			template<> struct has_fixed_size<T> {											\
				static const bool value = true;												\
			};																				\
			template<> struct fixed_size_helper<T, T> {										\
				static const std::size_t value = 1;											\
			};																				\
			template<> struct base_type<T> {												\
				typedef T type;																\
			};
		#define MOCASITO_IO_SET_SCALAR_FIXED_CV(T, N)										\
			MOCASITO_IO_SET_SCALAR_FIXED(T, N)												\
			MOCASITO_IO_SET_SCALAR_FIXED(T const, N)										\
			MOCASITO_IO_SET_SCALAR_FIXED(T volatile, N)										\
			MOCASITO_IO_SET_SCALAR_FIXED(T const volatile, N)								\
			MOCASITO_IO_GET_DATA(T, v, &v)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(bool, 1)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(char, 2)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(signed char, 3)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(unsigned char, 4)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(short, 5)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(unsigned short, 6)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(int, 7)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(unsigned, 8)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(long, 9)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(unsigned long, 10)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(long long, 11)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(unsigned long long, 12)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(float, 13)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(double, 14)
		MOCASITO_IO_SET_SCALAR_FIXED_CV(long double, 15)
		// complex types
		#define MOCASITO_IO_SET_FIXED(P, T, C, N)											\
			template<P> struct has_fixed_size<T> {											\
				static const bool value = has_fixed_size<C>::value;							\
			};																				\
			template<P> struct value_traits<T> {											\
				typedef T type;																\
			};																				\
			template<P> struct fixed_size_helper<T, T> {									\
				static const std::size_t value = N * fixed_size<C>::value;					\
			};																				\
			template<P> struct base_type<T> {												\
				typedef typename base_type<C>::type type;									\
			};
		#define MOCASITO_IO_SET_FIXED_CV(P, T, C, N)										\
			MOCASITO_IO_SET_FIXED(P, T, C, N)												\
			MOCASITO_IO_SET_FIXED(P, T const, C, N)											\
			MOCASITO_IO_SET_FIXED(P, T volatile, C, N)										\
			MOCASITO_IO_SET_FIXED(P, T const volatile, C, N)
		MOCASITO_IO_SET_FIXED_CV(typename T, std::complex<T>, T, 2)
		// C array
		template <typename T, typename Base> struct fixed_size_helper<T, Base[]> {
			static const size_t value = sizeof(T)/sizeof(Base);
		};
		template<typename T> struct has_fixed_size<T[]> {
			static const bool value = has_fixed_size<T>::value;
		};
		template<typename T> struct base_type<T[]> {
			typedef typename base_type<T>::type type;
		};
		template<typename T> struct value_traits<T[]> {
			typedef T type;
		};
		template<typename T> struct value_traits<T const []> {
			typedef T const type;
		};
		// std::string
		template<> struct value_traits<std::string> {
			typedef char type;
		};
		MOCASITO_IO_GET_DATA(std::string, v, &v[0])
	}
}
#endif

