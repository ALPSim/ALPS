/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_PARAMVALUE_READER_HPP
#define ALPS_NGS_DETAIL_PARAMVALUE_READER_HPP

#include <alps/ngs/cast.hpp>
#include <alps/ngs/config.hpp>

#include <boost/variant.hpp>

namespace alps {
    namespace detail {

        template<typename T> struct paramvalue_reader_visitor {
            
            template <typename U> void operator()(U const & data) {
                value = cast<T>(data);
            }
            
            template <typename U> void operator()(U * const, std::vector<std::size_t>) {
                throw std::runtime_error(std::string("cannot cast from std::vector<") + typeid(U).name() + "> to " + typeid(T).name() + ALPS_STACKTRACE);
            }

            T value;
        };

        template<typename T> struct paramvalue_reader_visitor<std::vector<T> > {
            
            template <typename U> void operator()(U const & data) {
                value.push_back(cast<T>(data));
            }

            template <typename U> void operator()(U * const ptr, std::vector<std::size_t> size) {
                if (size.size() != 1)
                    throw std::invalid_argument("only 1 D array are supported in alps::params" + ALPS_STACKTRACE);
                else if (size[0] != 0)
                    for (U const * it = ptr; it != ptr + size[0]; ++it)
                        (*this)(*it);
            }

            std::vector<T> value;
        };

        template<> struct paramvalue_reader_visitor<std::string> {
            
            template <typename U> void operator()(U const & data) {
                value = cast<std::string>(data);
            }
            
            template <typename U> void operator()(U * const ptr, std::vector<std::size_t> size) {
                if (size.size() != 1)
                    throw std::invalid_argument("only 1 D array are supported in alps::params" + ALPS_STACKTRACE);
                else if (size[0] != 0)
                    for (U const * it = ptr; it != ptr + size[0]; ++it)
                        value += (it == ptr ? "," : "") + cast<std::string>(*it);
            }

            std::string value;
        };

        template<typename T> struct paramvalue_reader 
            : public boost::static_visitor<> 
        {
            public:

                template <typename U> void operator()(U const & v) const {
                    visitor(v);
                }
                
                template <typename U> void operator()(std::vector<U> const & v) const {
                    visitor(v.data(), std::vector<std::size_t>(1, v.size()));
                }

                // std::vector<bool> stores proxy bits rather than contiguous
                // bool objects and therefore has no usable data() pointer.
                // Materialise byte values for the existing conversion visitor;
                // scalar targets still reject vector input, while vector targets
                // convert each byte to their requested element type.
                void operator()(std::vector<bool> const & v) const {
                    std::vector<unsigned char> contiguous(v.begin(), v.end());
                    visitor(contiguous.data(),
                            std::vector<std::size_t>(1, contiguous.size()));
                }

                T const & get_value() {
                    return visitor.value;
                }

            private:

                mutable paramvalue_reader_visitor<T> visitor;
        };

    }
}

#endif
