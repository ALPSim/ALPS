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

#include <iostream>
#include <alps/ngs.hpp>
#include <alps/ngs/numeric/array.hpp>
#include <alps/multi_array.hpp>

using namespace std;
using namespace alps::alea;
using namespace alps::ngs::numeric;

int main()
{
    #ifdef ALPS_NGS_USE_NEW_ALEA
        accumulator_set set;
        set << alps::ngs::RealObservable("Energy");
        set << alps::ngs::RealVectorObservable("Vel");
        set << alps::ngs::SimpleRealObservable("simple");
        
        set["Energy"] << 1.;
        set["Energy"] << 2.;
        set["Energy"] << 3.;
        set["Energy"] << 3.;
        
        std::vector<double> v1(3,1);
        std::vector<double> v2(3,2);
        std::vector<double> v3(3,3);
        std::vector<double> v4(3,3);
        
        set["Vel"] << v1;
        set["Vel"] << v2;
        set["Vel"] << v3;
        set["Vel"] << v4;
        
        cout << "mean: " << set["Energy"].get<double>().mean() << endl;
        cout << "mean: " << alps::short_print(set["Vel"].get<std::vector<double> >().mean()) << endl;
        cout << "error: " << alps::short_print(set["Vel"].get<std::vector<double> >().error()) << endl;
        
        alps::mcresult res(set["Energy"]);
        alps::mcresult res2(set["Vel"]);
        
        //~ alps::alea::accumulator<int> accc;
        //~ alps::alea::detail::accumulator_wrapper wa(accc);
        //~ alps::mcresult res3(wa);
        
        set["Vel"].get<std::vector<double> >();
        
        std::cout << res << std::endl;
        std::cout<< res2 << std::endl;
        
        //------------------- test boost::array -------------------
        
        boost::array<double, 3> a;
        
        for(int i = 0; i < 3; ++i)
        {
            a[i] = i;
        }
        alps::alea::accumulator<boost::array<double, 3>, features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                //~ , tag::converged
                                //~ , tag::tau
                                , tag::histogram
                                > > accc;
        
        accc << a;                        
        
        alps::alea::detail::accumulator_wrapper wa(accc);
        
        set << alps::alea::make_accumulator(accc, "bar");
        
        //------------------- test reset -------------------
        set["bar"] << a;
        std::cout << "Mean: " << alps::short_print(set["bar"].get<boost::array<double, 3> >().mean()) << std::endl;
        set.reset();
        set["bar"] << a+a;
        std::cout << "Mean: " << alps::short_print(set["bar"].get<boost::array<double, 3> >().mean()) << std::endl;
        
        //------------------- test array ops -------------------
        a+=sqrt(a+a-a*a/2);
        
        boost::array<double, 3> b = boost::array<double, 3>();
        for(int i = 0; i < 3; ++i)
            std::cout << b[i] << std::endl;

        b[0]=1;
        b[1]=1;
        b[2]=1;
        
        accc << b;
        
        for(int i = 0; i < 3; ++i)
        {
            //~ std::cout << a[i] << std::endl;
        }
        //~ std::cout << alps::short_print(accc.mean()) << std::endl;
        
        //------------------- test boost::multi-array -------------------
        alps::alea::accumulator<alps::multi_array<double, 3>, features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                //~ , tag::converged
                                //~ , tag::tau
                                , tag::histogram
                                > > accma;
        
        alps::multi_array<double, 3> ma(2,2,2);

        for(uint i = 0; i < 2; ++i)
        {
            for(uint j = 0; j < 2; ++j)
            {
                for(uint k = 0; k < 2; ++k)
                {
                    ma[i][j][k] = i*4+j*2+k;
                }
            }
        }
        
        accma << ma;
        
        std::cout << alps::short_print(ma) << std::endl;
        
    #else
        std::cout << "ALPS_NGS_USE_NEW_ALEA is OFF" << std::endl;
    #endif
}
