/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_ALEA_RESULT_TYPE_WRAPPER_HPP
#define ALPS_NGS_ALEA_RESULT_TYPE_WRAPPER_HPP

#include <alps/ngs/alea/features.hpp>
#include <alps/ngs/alea/wrapper/base_wrapper.hpp>

namespace alps {
    namespace accumulator {
        namespace detail {

            template <typename ValueType> class result_type_result_wrapper: public base_result_wrapper {
                public:
                    typedef ValueType value_type;
                    virtual ~result_type_result_wrapper() {}
                    virtual typename mean_type<value_type>::type mean() const = 0;
                    virtual bool has_mean() const = 0;
                    virtual typename error_type<value_type>::type error() const = 0;
                    virtual bool has_error() const = 0;
                    virtual typename tau_type<value_type>::type tau() const = 0;
                    virtual bool has_tau() const = 0;
            };            

            template <typename ValueType> class result_type_accumulator_wrapper: public base_accumulator_wrapper {
                public:
                    typedef ValueType value_type;
                    virtual ~result_type_accumulator_wrapper() {}
                    virtual typename mean_type<value_type>::type mean() const = 0;
                    virtual bool has_mean() const = 0;
                    virtual typename error_type<value_type>::type error() const = 0;
                    virtual bool has_error() const = 0;
                    virtual typename fixed_size_binning_type<value_type>::type fixed_size_binning() const = 0;
                    virtual bool has_fixed_size_binning() const = 0;
                    virtual typename max_num_binning_type<value_type>::type max_num_binning() const = 0;
                    virtual bool has_max_num_binning() const = 0;
                    virtual typename log_binning_type<value_type>::type log_binning() const = 0;
                    virtual bool has_log_binning() const = 0;
                    virtual typename autocorrelation_type<value_type>::type autocorrelation() const = 0;
                    virtual bool has_autocorrelation() const = 0;
                    virtual typename converged_type<value_type>::type converged() const = 0;
                    virtual bool has_converged() const = 0;
                    virtual typename tau_type<value_type>::type tau() const = 0;
                    virtual bool has_tau() const = 0;
                    virtual typename weight_type<value_type>::type weight() const = 0;
                    virtual bool has_weight() const = 0;
                    virtual typename histogram_type<value_type>::type histogram() const = 0;
                    virtual bool has_histogram() const = 0;
            };
        }
    }
}
#endif
