/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef UTILS_DMRG_RANDOM_HPP
#define UTILS_DMRG_RANDOM_HPP

#include <boost/random.hpp>

struct dmrg_random {
    typedef double value_type;
    typedef boost::mt19937 engine_t;
    typedef boost::uniform_real<value_type> uniform_dist_t;
    typedef boost::normal_distribution<value_type> normal_dist_t;
    typedef boost::poisson_distribution<value_type> poisson_dist_t;
    
    static engine_t engine;

    // Uniform distribution
    static inline value_type uniform (value_type min, value_type max) {
        uniform_dist_t dist(min, max);
        return dist(engine);
    }
    
    static inline value_type uniform () {
        return uniform(0, 1);
    }

    
    // Normal distribution
    static inline value_type normal (value_type mean, value_type sigma) {
        normal_dist_t dist(mean, sigma);
        return dist(engine);
    }
    
    static inline value_type normal () {
        return normal(0, 1);
    }

    
    // Poisson distribution
    /*
    static inline value_type poisson (value_type mean) {
        poisson_dist_t dist(mean);
        return dist(engine);
    }
    
    static inline value_type poisson () {
        return poisson(1);
    }
*/
    
};

#endif
