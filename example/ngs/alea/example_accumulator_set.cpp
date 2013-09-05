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
#include <alps/multi_array.hpp>
#include <alps/utility/resize.hpp>
#include <alps/hdf5/array.hpp>
#include <alps/ngs/numeric/array.hpp>

#include <iostream>

using namespace std;
using namespace alps::accumulator;
using namespace alps::ngs::numeric;

int main() {
    #ifdef ALPS_NGS_USE_NEW_ALEA
        
        //=================== test resize_same_as ===================
        std::cout << "resize_same_as test" << std::endl;
        //~ typedef alps::multi_array<double, 3> exp_type;
        //~ exp_type a(2,2,2);
        //~ typedef std::vector<double> exp_type;        
        //~ exp_type a(2, 3);
        typedef boost::array<double, 3> exp_type;
        exp_type a;
        
        exp_type b;
        
        alps::resize_same_as(b, a);
        
        std::cout << alps::short_print(a) << std::endl;
        std::cout << alps::short_print(b) << std::endl;
        
        //=================== test deep copy ===================
        std::cout << "deep_copy test" << std::endl;
        
        accumulator<double> original;
        original << 1.;
        alps::accumulator::detail::accumulator_wrapper copy1(original);
        original << 1.;
        copy1 << 3.;
        accumulator_set copy2;
        copy2 << make_accumulator("copy", original);
        original << 1.;
        copy2["copy"] << 7.;
        
        std::cout << mean(original) << std::endl;
        std::cout << mean(copy1.get<double>()) << std::endl;
        std::cout << mean(copy2["copy"].get<double>()) << std::endl;
        
        //=================== the supported Observables ===================
        std::cout << "set test" << std::endl;
        accumulator_set set;
        set << alps::ngs::RealObservable("Energy");
        set << alps::ngs::RealVectorObservable("Vel");
        set << alps::ngs::SimpleRealObservable("simple");
        set << alps::ngs::SimpleRealVectorObservable("simpleV");
        
        //------------------- stream in values (must be of value_type) -------------------
        set["Energy"] << 1.; // ... << 1; wouldn't work
        set["Energy"] << 2.;
        set["Energy"] << 3.;
        set["Energy"] << 3.;
        
        set["simple"] << 2.;
        set["simple"] << 3.;
        
        set["Vel"] << std::vector<double>(3,1);
        set["Vel"] << std::vector<double>(3,2);
        set["Vel"] << std::vector<double>(3,3);
        set["Vel"] << std::vector<double>(3,4);
        
        set["simpleV"] << std::vector<double>(3,1);
        set["simpleV"] << std::vector<double>(3,2);
        
        //------------------- get some data via the result_type_accumulator_wrapper -------------------
        //------------------- manual type-infusion -------------------
        cout << "mean: " << set["Energy"].get<double>().mean() << endl;
        cout << "mean: " << alps::short_print(set["Vel"].get<std::vector<double> >().mean()) << endl;
        cout << "error: " << alps::short_print(set["Vel"].get<std::vector<double> >().error()) << endl;
        
        //=================== convert to mcresult ===================
        //------------------- automatic type-infusion (try...catch) -------------------
        alps::mcresult res(set["Energy"]);
        alps::mcresult res_v(set["Vel"]);
        alps::mcresult res_s(set["simple"]);
        alps::mcresult res_s_v(set["simpleV"]);
        
        //------------------- any accumulator can be inserted in the set with make_accumulator -------------------
        alps::accumulator::accumulator<double, features<tag::mean> > acc;
        //~ detail::accumulator_wrapper wa0(acc);
        //~ wa0 << 1.;
        
        set << make_accumulator("some_wrapper", acc);
        
        
        set["some_wrapper"] << 2.;
        set["some_wrapper"] << 3.;
        //------anything in a wrapper can be casted to mcresult (only some value_type work (double, ...)) ---------
        //~ alps::mcresult res_wrapper(wa0);
        alps::mcresult res_wrapper(set["some_wrapper"]);
        
        //------------------- control-print -------------------
        std::cout << "RealObservable:             " << res << std::endl;
        std::cout << "RealVectorObservable:       " << res_v << std::endl;
        std::cout << "SimpleRealObservable:       " << res_s << std::endl;
        std::cout << "SimpleRealVectorObservable: " << res_s_v << std::endl;
        std::cout << "Some wrapper:               " << res_wrapper << std::endl;
        
        //=================== test boost::array ===================
        boost::array<double, 3> boost_array;
        
        for(int i = 0; i < 3; ++i)
        {
            boost_array[i] = i;
        }
        alps::accumulator::accumulator<boost::array<double, 3>, features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                , tag::histogram
                                > > acc_boost_array;
        
        acc_boost_array << boost_array;
        
        
        set << make_accumulator("bar", acc_boost_array);
        
        //------------------- test reset -------------------
        set["bar"] << boost_array;
        std::cout << "Mean: " << alps::short_print(set["bar"].get<boost::array<double, 3> >().mean()) << std::endl;
        set.reset();
        set["bar"] << boost_array + boost_array;
        std::cout << "Mean: " << alps::short_print(set["bar"].get<boost::array<double, 3> >().mean()) << std::endl;
        
        
        //------------------- doesn't work yet -------------------
        //~ alps::mcresult res4(set["bar"]);
        //~ std::cout << res4 << std::endl;
        
        //------------------- test array ops -------------------
        //------------------- they fulfill the error_type_concept -------------------
        boost_array += sqrt(boost_array + boost_array - boost_array * boost_array / 2.);
        
        //=================== test alps::multi-array ===================
        alps::accumulator::accumulator<alps::multi_array<double, 3>, features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                , tag::histogram
                                > > acc_alps_ma;
        //TODO: multiarray not supported yet. See check_size in ngs/numeric/detail.hpp
        //~ alps::multi_array<double, 3> ma(2,2,2);
        //~ 
        //~ for(uint i = 0; i < 2; ++i)
        //~ {
            //~ for(uint j = 0; j < 2; ++j)
            //~ {
                //~ for(uint k = 0; k < 2; ++k)
                //~ {
                    //~ ma[i][j][k] = i*4+j*2+k;
                //~ }
            //~ }
        //~ }
        
        //~ acc_alps_ma << ma;
        
        //~ std::cout << alps::short_print(ma) << std::endl;
        
    #else
        std::cout << "ALPS_NGS_USE_NEW_ALEA is OFF" << std::endl;
    #endif
}
