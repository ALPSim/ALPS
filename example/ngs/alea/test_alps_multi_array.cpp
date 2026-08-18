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

#include <alps/ngs/short_print.hpp>
#include <alps/multi_array.hpp>

#include <iostream>
#include <vector>

#include <boost/integer.hpp>

using namespace std;

int main()
{
    alps::multi_array<double, 3> a(4,2,2);
    alps::multi_array<double, 2> b(4,4);
    alps::multi_array<double, 3> c(4,2,2);
    
    for(int i = 0; i < 4; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            for(int k = 0; k < 2; ++k)
            {
                a[i][j][k] = 4*i+2*j+k;
                c[i][j][k] = 4*(i%2)+2*j+k;
            }
        }
    }
    
    
    a+a;
    std::cout << sqrt(a) << std::endl;
    std::cout << "--" << std::endl;
    std::cout << a-c << std::endl;
    std::cout << "--" << std::endl;
    c = a;
    std::cout << alps::sqrt(a) << std::endl;
}
