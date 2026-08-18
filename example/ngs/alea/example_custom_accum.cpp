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


#include <alps/ngs.hpp>

#include "custom_accum.hpp"

#include <iostream>

using namespace std;
using namespace alps::accumulator;

int main()
{
    cout << "test custom_accum" << endl;
    cout << "-----------------" << endl;
    
    custom_accum a;
    
    add_value(a, 1);
    add_value(a, 2);
    add_value(a, 3);
    add_value(a, 4);
    add_value(a, 5);
    
    count(a);
    mean(a);
    
    cout << a << endl;
    
    detail::accumulator_wrapper m(a);
    
    m(6);
    m(7);
    m(8);
    
    // get infos via intermediate type m.get<int>(), where int is the value_type of custom_accum
    m.get<int>().count();
    m.get<int>().mean();
    
    // extract the accumulator (and get count/mean info)
    count(extract<custom_accum>(m));
    mean(m.extract<custom_accum>());
    
    cout << m;
}
