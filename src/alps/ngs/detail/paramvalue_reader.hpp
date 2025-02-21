/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_PARAMVALUE_READER_HPP
#define ALPS_NGS_DETAIL_PARAMVALUE_READER_HPP

#include <alps/ngs/cast.hpp>
#include <alps/ngs/config.hpp>

#include <boost/variant.hpp>

#if defined(ALPS_HAVE_PYTHON)
	#include <alps/ngs/detail/get_numpy_type.hpp>
	#include <alps/ngs/detail/extract_from_pyobject.hpp>

	#include <alps/ngs/boost_python.hpp>
#endif

namespace alps {
    namespace detail {

        template<typename T> struct paramvalue_reader_visitor {
            
            template <typename U> void operator()(U const & data) {
                value = cast<T>(data);
            }
            
            template <typename U> void operator()(U * const, std::vector<std::size_t>) {
                throw std::runtime_error(std::string("cannot cast from std::vector<") + typeid(U).name() + "> to " + typeid(T).name() + ALPS_STACKTRACE);
            }

            #if defined(ALPS_HAVE_PYTHON)
                void operator()(boost::python::list const &) {
                    throw std::runtime_error(std::string("cannot cast from boost::python::list ") + typeid(T).name() + ALPS_STACKTRACE);
                }

                void operator()(boost::python::dict const &) {
                    throw std::invalid_argument("python dict cannot be used in alps::params" + ALPS_STACKTRACE);
                }
            #endif

            T value;
        };

        template<typename T> struct paramvalue_reader_visitor<std::vector<T> > {
            
            template <typename U> void operator()(U const & data) {
                value.push_back(cast<T>(data));
            }

            template <typename U> void operator()(U * const ptr, std::vector<std::size_t> size) {
                if (size.size() != 1)
                    throw std::invalid_argument("only 1 D array are supported in alps::params" + ALPS_STACKTRACE);
                else
                    for (U const * it = ptr; it != ptr + size[0]; ++it)
                        (*this)(*it);
            }

            #if defined(ALPS_HAVE_PYTHON)
                void operator()(boost::python::list const & data) {
                    for(boost::python::ssize_t i = 0; i < boost::python::len(data); ++i) {
                        paramvalue_reader_visitor<T> scalar;
                        extract_from_pyobject(scalar, data[i]);
                        value.push_back(scalar.value);
                    }
                }

                void operator()(boost::python::dict const &) {
                    throw std::invalid_argument("python dict cannot be used in alps::params" + ALPS_STACKTRACE);
                }
            #endif

            std::vector<T> value;
        };

        template<> struct paramvalue_reader_visitor<std::string> {
            
            template <typename U> void operator()(U const & data) {
                value = cast<std::string>(data);
            }
            
            template <typename U> void operator()(U * const ptr, std::vector<std::size_t> size) {
                if (size.size() != 1)
                    throw std::invalid_argument("only 1 D array are supported in alps::params" + ALPS_STACKTRACE);
                else
                    for (U const * it = ptr; it != ptr + size[0]; ++it)
                        value += (it == ptr ? "," : "") + cast<std::string>(*it);
            }

            #if defined(ALPS_HAVE_PYTHON)
                void operator()(boost::python::list const & data) {
                    for(boost::python::ssize_t i = 0; i < boost::python::len(data); ++i)
                        value += (value.size() ? "," : "") + boost::python::call_method<std::string>(boost::python::object(data[i]).ptr(), "__str__");
                }

                void operator()(boost::python::dict const &) {
                    throw std::invalid_argument("python dict cannot be used in alps::params" + ALPS_STACKTRACE);
                }
            #endif

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
                    visitor(&v.front(), std::vector<std::size_t>(1, v.size()));
                }

                void operator()(T const & v) const {
                    visitor.value = v; 
                }

                #if defined(ALPS_HAVE_PYTHON)
                    void operator()(boost::python::object const & v) const {
                        extract_from_pyobject(visitor, v);
                    }
                #endif

                T const & get_value() {
                    return visitor.value;
                }

            private:

                mutable paramvalue_reader_visitor<T> visitor;
        };

        #if defined(ALPS_HAVE_PYTHON)
            template<> struct paramvalue_reader<boost::python::object>
                : public boost::static_visitor<> 
            {
                public:

                    template <typename U> void operator()(U const & v) const {
                        value = boost::python::object(v);
                    }

                    template <typename U> void operator()(std::vector<U> const & v) const {
                        npy_intp npsize = v.size();
                        value = boost::python::object(boost::python::handle<>(PyArray_SimpleNew(1, &npsize, detail::get_numpy_type(U()))));
                        PyArrayObject * ptr = (PyArrayObject *)value.ptr();
                        memcpy(PyArray_DATA(ptr), &v.front(), PyArray_ITEMSIZE(ptr) * PyArray_SIZE(ptr));
                    }

                    void operator()(std::vector<std::string> const & v) const {
                        value = boost::python::list(v);
                    }

                    void operator()(boost::python::object const & v) const {
                        value = v; 
                    }

                    boost::python::object const & get_value() {
                        return value;
                    }

                private:

                    mutable boost::python::object value;
            };
        #endif
    }
}

#endif
