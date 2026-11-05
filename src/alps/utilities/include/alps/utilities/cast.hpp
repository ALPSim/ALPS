/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_UTILITY_CAST_HPP
#define ALPS_UTILITY_CAST_HPP

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <type_traits>

#include <alps/ngs/config.hpp>
#include <alps/utilities/stacktrace.hpp>

#include <filesystem>

#include <string>
#include <complex>
#include <typeinfo>
#include <algorithm>
#include <stdexcept>
#include <cstdio>

namespace alps {

    struct bad_cast : public std::runtime_error {
        bad_cast(std::string const & arg) : std::runtime_error(arg) {}
    };

    template<typename U, typename T> inline U cast(T const &);

    namespace detail {
    
        template<typename U, typename T> struct is_cast {
            static T t;
            static char check(U);
            static double check(...);
            enum { value = sizeof(check(t)) / sizeof(char) };
        };

        template<
            typename U, typename T, typename X
        > inline U cast_generic(T arg, X) {
            throw bad_cast(
                  std::string("cannot cast from ") 
                + typeid(T).name() 
                + " to " 
                + typeid(U).name() + ALPS_STACKTRACE
            );
            return U();
        }

        template<typename U, typename T> inline U cast_generic(
            T arg, std::integral_constant<int, 1> const&
        ) {
            return arg;
        }

    }

    template<typename U, typename T> struct cast_hook {
        static inline U apply(T arg) {
            return detail::cast_generic<U, T>(
                arg, std::integral_constant<int, detail::is_cast<U, T>::value>()
            );
        }
    };

    namespace detail {
        // The text must be one number, apart from surrounding whitespace;
        // "6x", "2*L" or "" are not silently read as their numeric prefix.
        inline bool only_space(char const * rest) {
            while (*rest && std::isspace(static_cast<unsigned char>(*rest))) ++rest;
            return !*rest;
        }

        template <typename T> T checked_real_string(std::string const & text) {
            char const * begin = text.c_str();
            char * end = nullptr;
            T value;
            if constexpr (std::is_same<T, float>::value) value = std::strtof(begin, &end);
            else if constexpr (std::is_same<T, double>::value) value = std::strtod(begin, &end);
            else value = std::strtold(begin, &end);
            if (end == begin || !only_space(end))
                throw std::invalid_argument("parameter is not a number: \"" + text + "\"");
            return value;
        }

        template <typename T> T checked_integer_string(std::string const & text) {
            // Keep the historical decimal parsing (e.g. "4.0" -> 4), but
            // never let scanf silently wrap an out-of-range parameter.
            char const * begin = text.c_str();
            char * end = nullptr;
            errno = 0;
            if constexpr (std::is_signed<T>::value) {
                long long value = std::strtoll(begin, &end, 10);
                if (end != begin && !only_space(end)) {
                    long double real = checked_real_string<long double>(text);
                    if (!(real > (std::numeric_limits<long long>::min)() - 1.L
                          && real < (std::numeric_limits<long long>::max)() + 1.L))
                        throw std::out_of_range("integer parameter out of range: " + text);
                    value = static_cast<long long>(real);
                } else if (end == begin)
                    throw std::invalid_argument("parameter is not an integer: \"" + text + "\"");
                if (errno == ERANGE || value < (std::numeric_limits<T>::min)()
                                    || value > (std::numeric_limits<T>::max)())
                    throw std::out_of_range("integer parameter out of range: " + text);
                return static_cast<T>(value);
            } else {
                auto first = text.find_first_not_of(" \t\r\n\f\v");
                if (first != std::string::npos && text[first] == '-')
                    throw std::out_of_range("integer parameter out of range: " + text);
                unsigned long long value = std::strtoull(begin, &end, 10);
                if (end != begin && !only_space(end)) {
                    long double real = checked_real_string<long double>(text);
                    if (!(real < (std::numeric_limits<unsigned long long>::max)() + 1.L))
                        throw std::out_of_range("integer parameter out of range: " + text);
                    value = static_cast<unsigned long long>(real);
                } else if (end == begin)
                    throw std::invalid_argument("parameter is not an integer: \"" + text + "\"");
                if (errno == ERANGE || value > (std::numeric_limits<T>::max)())
                    throw std::out_of_range("integer parameter out of range: " + text);
                return static_cast<T>(value);
            }
        }

        inline bool checked_bool_string(std::string const & text) {
            std::string word;
            for (char c : text)
                if (!std::isspace(static_cast<unsigned char>(c)))
                    word += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (word == "true") return true;
            if (word == "false") return false;
            return checked_integer_string<short>(text) != 0;
        }
    }

    #define ALPS_CAST_STRING(T, p, c)                                            \
        template<> struct cast_hook<std::string, T > {                                \
            static inline std::string apply( T arg) {                                \
                char buffer[255];                                                    \
                if (snprintf(buffer, 255, "%" p "" c, arg) < 0)                       \
                    throw std::runtime_error(                                        \
                        "error casting from " #T " to string" + ALPS_STACKTRACE        \
                    );                                                                \
                return buffer;                                                        \
            }                                                                        \
        };                                                                            \
        template<> struct cast_hook< T, std::string> {                                \
            static inline T apply(std::string arg) {                                \
                if constexpr (std::is_integral<T>::value)                           \
                    return detail::checked_integer_string<T>(arg);                 \
                else                                                                \
                    return detail::checked_real_string<T>(arg);                    \
            }                                                                        \
        };
    ALPS_CAST_STRING(short, "", "hd")
    ALPS_CAST_STRING(int, "", "d")
    ALPS_CAST_STRING(long, "", "ld")
    ALPS_CAST_STRING(unsigned short, "", "hu")
    ALPS_CAST_STRING(unsigned int, "", "u")
    ALPS_CAST_STRING(unsigned long, "", "lu")
    ALPS_CAST_STRING(float, ".8", "e")
    ALPS_CAST_STRING(double, ".16", "le")
    ALPS_CAST_STRING(long double, ".32", "Le")
    ALPS_CAST_STRING(long long, "", "lld")
    ALPS_CAST_STRING(unsigned long long, "", "llu")
    #undef ALPS_CAST_STRING

    #define ALPS_CAST_STRING_CHAR(T, U)                                            \
        template<> struct cast_hook<std::string, T > {                                \
            static inline std::string apply( T arg) {                                \
                return cast_hook<std::string, U>::apply(arg);                        \
            }                                                                        \
        };                                                                            \
        template<> struct cast_hook<T, std::string> {                                \
            static inline T apply(std::string arg) {                                \
                return cast_hook< U , std::string>::apply(arg);                        \
            }                                                                        \
        };
    template<> struct cast_hook<std::string, bool> {
        static inline std::string apply(bool arg) {
            return cast_hook<std::string, short>::apply(arg);
        }
    };
    template<> struct cast_hook<bool, std::string> {
        static inline bool apply(std::string arg) {
            return detail::checked_bool_string(arg);
        }
    };
    ALPS_CAST_STRING_CHAR(char, short)
    ALPS_CAST_STRING_CHAR(signed char, short)
    ALPS_CAST_STRING_CHAR(unsigned char, unsigned short)
    #undef ALPS_CAST_STRING_CHAR

    template<typename U, typename T> struct cast_hook<U, std::complex<T> > {
        static inline U apply(std::complex<T> const & arg) {
            return static_cast<U>(arg.real());
        }
    };

    template<typename U, typename T> struct cast_hook<std::complex<U>, T> {
        static inline std::complex<U> apply(T const & arg) {
            return cast<U>(arg);
        }
    };

    template<typename U, typename T> struct cast_hook<std::complex<U>, std::complex<T> > {
        static inline std::complex<U> apply(std::complex<T> const & arg) {
            return std::complex<U>(arg.real(), arg.imag());
        }
    };

    template<typename T> struct cast_hook<std::string, std::complex<T> > {
        static inline std::string apply(std::complex<T> const & arg) {
            return cast<std::string>(arg.real()) + "+" + cast<std::string>(arg.imag()) + "i";
        }
    };

    // TODO: also parse a+bi
    template<typename T> struct cast_hook<std::complex<T>, std::string> {
        static inline std::complex<T> apply(std::string const & arg) {
            return cast<T>(arg);
        }
    };
    
    template<typename T> struct cast_hook<std::filesystem::path, T> {
        static inline std::filesystem::path apply(T const & arg) {
            return std::filesystem::path(cast<std::string>(arg));
        }
    };

    template<typename T> struct cast_hook<std::filesystem::path, std::complex<T> > {
        static inline std::filesystem::path apply(std::complex<T> const & arg) {
            return std::filesystem::path(cast<std::string>(arg));
        }
    };

    template<typename U, typename T> inline U cast(T const & arg) {
        return cast_hook<U, T>::apply(arg);
    }

    template<typename U, typename T> inline void cast(
        U const * src, U const * end, T * dest
    ) {
        for (U const * it = src; it != end; ++it)
            dest[it - src] = cast<T>(*it);
    }
}

#endif
