/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Lukas Gamper <gamperl@gmail.ch>                           *
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

#ifndef ALPS_NGS_ALEA_RESULT_WRAPPER_HPP
#define ALPS_NGS_ALEA_RESULT_WRAPPER_HPP

#include <alps/ngs/stacktrace.hpp>
#include <alps/ngs/alea/result.hpp>
#include <alps/ngs/alea/wrapper/base_wrapper.hpp>
#include <alps/ngs/alea/wrapper/derived_wrapper.hpp>
#include <alps/ngs/alea/wrapper/result_type_wrapper.hpp>
 
#ifdef ALPS_HAVE_MPI
    #include <alps/ngs/boost_mpi.hpp>
#endif

#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>

#include <typeinfo> //used in add_value
#include <stdexcept>

namespace alps{
    namespace accumulator {
        namespace detail {

            // class that holds the base_result_wrapper pointer
            class result_wrapper {
                public:
                    result_wrapper(boost::shared_ptr<base_result_wrapper> const & arg)
                        : base_(arg)
                    {}

                    result_wrapper(result_wrapper const & arg)
                        : base_(arg.base_->clone())
                    {}

                    template<typename T> result_type_result_wrapper<T> & get() const {
                        return base_->get<T>();
                    }

                    friend std::ostream& operator<<(std::ostream &out, result_wrapper const & wrapper);

                    template <typename T> T & extract() const {
                        return (dynamic_cast<derived_result_wrapper<T>& >(*base_)).accum_;
                    }

                    boost::uint64_t count() const {
                        return base_->count();
                    }

                private:
                    boost::shared_ptr<base_result_wrapper> base_;
            };

            inline std::ostream & operator<<(std::ostream & out, result_wrapper const & m) {
                m.base_->print(out);
                return out;
            }
        }

        template <typename Result> Result & extract(detail::result_wrapper & m) {
            return m.extract<Result>();
        }
    }
}
#endif
