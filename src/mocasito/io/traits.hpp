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
		#define MOCASITO_IO_CV_CALLBACK_1(callback, T, arg0)								\
			callback(T, arg0)																\
			callback(T const, arg0)															\
			callback(T volatile, arg0)														\
			callback(T const volatile, arg0)
		template<typename T = void> struct type_traits {
			typedef int type;
			static const int value = -1;
		};
		template<typename T = void> struct pointer_traits {
			typedef typename pointer_traits<typename T::value_type>::type type;
		};
		template<typename T> struct pointer_traits<T const> {
			typedef typename pointer_traits<typename T::value_type const>::type type;
		};
		template<typename T> struct base_type {
			typedef typename base_type<typename T::value_type>::type type;
		};
		template<typename T> struct base_type<T const> {
			typedef typename base_type<typename T::value_type>::type const type;
		};
		template <typename T> struct fixed_size {
			static std::size_t const value = 0;
		};
		template<typename T> struct get_data {
			static typename pointer_traits<T>::type apply(T & v) {
				return get_data<typename T::value_type>::apply(*v.begin());
			}
		};
		template<typename T> struct get_data<T const> {
			static typename pointer_traits<T const>::type apply(T const & v) {
				return get_data<typename T::value_type const>::apply(*v.begin());
			}
		};
		template<typename T> struct get_size {
			static std::size_t apply(T & v) {
				return v.size() * fixed_size<typename T::value_type>::value;
			}
		};
		template<typename T> struct set_size {
			static void apply(T & v, std::size_t s) {
				v.resize(s / fixed_size<typename T::value_type>::value);
			}
		};
		// SCALAR
		#define MOCASITO_IO_SET_SCALAR(T, N)												\
			template<> struct type_traits<T> {												\
				typedef int type;															\
				static int const value = N;													\
			};																				\
			template<> struct pointer_traits<T> {											\
				typedef T * type;															\
			};																				\
			template<> struct base_type<T> {												\
				typedef T type;																\
			};																				\
			template <> struct fixed_size<T> {												\
				static std::size_t const value = 1;											\
			};																				\
			template<> struct get_data<T> {													\
				static pointer_traits<T>::type apply(T & v) {								\
					return &v;																\
				}																			\
			};																				\
			template<> struct get_size<T> {													\
				static std::size_t apply(T &) {											\
					return fixed_size<base_type<T>::type>::value;							\
				}																			\
			};																				\
			template<> struct set_size<T> {													\
				static void apply(T, std::size_t) {}										\
			};
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, bool, 1)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, char, 2)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, signed char, 3)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, unsigned char, 4)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, short, 5)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, unsigned short, 6)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, int, 7)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, unsigned, 8)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, long, 9)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, unsigned long, 10)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, long long, 11)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, unsigned long long, 12)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, float, 13)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, double, 14)
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_SCALAR, long double, 15)
		// COMPLEX
		template<typename T> struct pointer_traits<std::complex<T> > {
			typedef typename pointer_traits<T>::type type;
		};
		template<typename T> struct pointer_traits<std::complex<T> const> {
			typedef typename pointer_traits<const T>::type type;
		};
		#define MOCASITO_IO_SET_COMPLEX(C, unused)											\
			template<typename T> struct base_type<C> {										\
				typedef typename base_type<T>::type type;									\
			};																				\
			template<typename T> struct type_traits<C> {									\
				static int const value = type_traits<typename base_type<C>::type>::value;	\
			};																				\
			template <typename T> struct fixed_size<C> {									\
				static std::size_t const value = 2 * fixed_size<T>::value;					\
			};																				\
			template<typename T> struct get_data<C> {										\
				static typename pointer_traits<C>::type apply(C & v) {						\
					return &v.real();														\
				}																			\
			};																				\
			template<typename T> struct get_size<C> {										\
				static std::size_t apply(C) {												\
					return fixed_size<C>::value;											\
				}																			\
			};																				\
			template<typename T> struct set_size<C> {										\
				static void apply(T, std::size_t) {}										\
			};
		MOCASITO_IO_CV_CALLBACK_1(MOCASITO_IO_SET_COMPLEX, std::complex<T>, ~)
		#undef MOCASITO_IO_SET_COMPLEX
		// ARRAY
		template<typename T, std::size_t N> struct pointer_traits<T[N]> {
			typedef typename pointer_traits<T>::type type;
		};
		template<typename T, std::size_t N> struct pointer_traits<T const [N]> {
			typedef typename pointer_traits<const T>::type type;
		};
		template<typename T, std::size_t N> struct base_type<T[N]> {
			typedef typename base_type<T>::type type;
		};
		template <typename T, std::size_t N> struct fixed_size<T[N]> {
			static std::size_t const value = N * fixed_size<T>::value;
		};
		template<typename T, std::size_t N> struct get_data<T[N]> {
			static typename pointer_traits<T>::type apply(T * v) {
				return get_data<T>::apply(*v);
			}
		};
		template<typename T, std::size_t N> struct get_data<T const [N]> {
			static typename pointer_traits<T const>::type apply(T const * v) {
				return get_data<T const>::apply(*v);
			}
		};
		template<typename T, std::size_t N> struct get_size<T[N]> {
			static std::size_t apply(T[N]) {
				return fixed_size<T[N]>::value;
			}
		};
		template<typename T, std::size_t N> struct set_size<T[N]> {
			static void apply(T[N], std::size_t) {}
		};
		#undef MOCASITO_IO_CV_CALLBACK_1
	}
}
#endif

