/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <alps/ngs.hpp>

//these two flags will create the int main() together with unit_test.hpp
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
  
BOOST_AUTO_TEST_CASE(test_ctor_in_modular_accum)
{
    typedef alps::alea::accumulator<int
                                    , alps::alea::FixSizeBinning
                                    , alps::alea::MaxNumberBinning
                                    , alps::alea::Autocorrelation
                                    > accum;
    accum acc(alps::alea::bin_number = 10, alps::alea::bin_size = 10);
    
    acc << 1;
    acc << 2;
    acc << 3;
    acc << 4;
    acc << 5;
    
    accum acc2(acc);
        
    BOOST_REQUIRE( alps::alea::count(acc2) == acc.count());
    BOOST_REQUIRE( alps::alea::mean(acc2) == acc.mean());
    BOOST_REQUIRE( alps::alea::error(acc2) == acc.error());
    BOOST_REQUIRE( alps::alea::fix_size_bin(acc2) == acc.fix_size_bin());
    BOOST_REQUIRE( alps::alea::max_num_bin(acc2) == acc.max_num_bin());
    BOOST_REQUIRE( alps::alea::autocorr(acc2) == acc.autocorr());
    
}
