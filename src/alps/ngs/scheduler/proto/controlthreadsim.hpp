/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#if !defined(ALPS_NGS_SCHEDULER_CONTROLTHREADSIM_NG_HPP) && !defined(ALPS_NGS_SINGLE_THREAD)
#define ALPS_NGS_SCHEDULER_CONTROLTHREADSIM_NG_HPP

#include <alps/ngs/api.hpp>

#include <boost/thread.hpp>

namespace alps {

    template<typename Impl> class controlthreadsim_ng : public Impl {
        public:
            controlthreadsim_ng(typename alps::parameters_type<Impl>::type const & p, std::size_t seed_offset = 0)
                : Impl(p, seed_offset)
                , m_status(Impl::initialized)
            {
                // TODO: this is ugly, but virtual functions do not work in constructor of base class
                Impl::data_mutex = mutex(new native_lockable());
                Impl::result_mutex = mutex(new native_lockable());
            }

            double fraction_completed() const {
                typename Impl::lock_guard data_lock(Impl::get_data_lock());
                return Impl::fraction_completed();
            }

            typename Impl::status_type status() const {
                return m_status;
            }
        
        protected:

            template<typename T> class atomic {
                public:

                    atomic() {}
                    atomic(T const & v): value(v) {}
                    atomic(atomic<T> const & v): value(v.value) {}

                    atomic<T> & operator=(T const & v) {
                        boost::lock_guard<boost::mutex> lock(atomic_mutex);
                        value = v;
                        return *this;
                    }

                    operator T() const {
                        boost::lock_guard<boost::mutex> lock(atomic_mutex);
                        return value;
                    }

                private:

                    T volatile value;
                    boost::mutex mutable atomic_mutex;
            };

            void on_unlock() {}

            void set_status(typename Impl::status_type status) {
                m_status = status;
            }

        private:

            atomic<typename Impl::status_type> m_status;
    };

}

#endif
