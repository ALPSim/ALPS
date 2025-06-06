/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
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

#define BOOST_TEST_MODULE alps::ngs::accumulator

#include <alps/ngs.hpp>

#ifndef ALPS_LINK_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_fixed_size_bin_in_modular_accum)
{
    alps::accumulator::accumulator<int, alps::accumulator::features<alps::accumulator::tag::fixed_size_binning> > acci; //Default 128
    
    acci << 2;
    acci << 6;
    
    BOOST_REQUIRE( alps::accumulator::fixed_size_binning(acci).bin_size() == 128);
    
    
    alps::accumulator::accumulator<double, alps::accumulator::features<alps::accumulator::tag::fixed_size_binning> > accd(alps::accumulator::bin_size = 10);

    for(int i = 0; i < 100; ++i)
    {
        accd << 0.1*i;
    }
    
        
        
    BOOST_REQUIRE( alps::accumulator::fixed_size_binning(accd).bin_size() == 10);
    alps::accumulator::fixed_size_binning(accd);
    accd.fixed_size_binning();
}
