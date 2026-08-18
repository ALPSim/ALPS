/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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
