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

#include <alps/ngs/numeric/vector.hpp>
#include <alps/ngs/numeric/array.hpp>
#include <alps/multi_array.hpp>
#include <alps/ngs/short_print.hpp>
//~ #include <alps/ngs.hpp>

#include <boost/integer.hpp>

#include <iostream>
#include <vector>

using namespace std;
//~ using namespace alps::accumulator;
using namespace alps::ngs::numeric;

int main()
{
    std::vector<double> a(100,4);
    std::vector<double> a1(100,5);
    std::vector<double> a3(100,7);
    std::vector<double> b(100,2);
    std::cout << "a: " << alps::short_print(a) << std::endl;

    a+=a1;
    a/3.;
    std::cout << "a: " << alps::short_print(sqrt(a*b)) << std::endl;


    //~ accumulator<std::vector<double>, features<tag::mean, tag::fixed_size_binning> > c(bin_size = 1);
    //~ c << a;
    //~ std::cout << alps::short_print(c.mean()) << std::endl;
    //~ c << a1;
    //~ std::cout << alps::short_print(c.mean()) << std::endl;
    //~ c << a3;
    //~ std::cout << alps::short_print(c.mean()) << std::endl;
    //~ 
    //~ std::cout << alps::short_print(c.fixed_size_binning().bins()) << std::endl;
    //~ 
    //~ 
    //~ std::cout << c.count() << std::endl;
  
    boost::array<double, 3> ba;
    ba[0] = 3;
    ba[1] = 1;
    ba[2] = 2;
    ba += ba;
    std::cout << alps::short_print(sqrt(ba + ba)) << std::endl;
    
    alps::multi_array<double, 3> ma(boost::extents[3][2][2]);
    alps::multi_array<double, 3> mb;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            for(int k = 0; k < 2; ++k)
            {
                ma[i][j][k] = i*4+j*2+k;
            }
        }
    }
    
    mb+= ma;
    std::cout << alps::short_print(mb) << std::endl;
    
    
}
