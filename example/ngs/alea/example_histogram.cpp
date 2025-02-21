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

#include <alps/ngs.hpp>

#include <iostream>

using namespace std;
using namespace alps::accumulator;

int main()
{
    cout << "test int histo" << endl;
    cout << "--------------" << endl;
    
    /*                                                              updated: 18.12.2011
     * in the histogram_old the first template parameter stands for the 
     * value_type and the second for the weight_type
     * 
     * histo[value_type(1)] += weight_type(1)
     * adds weight_type(1) to the interval where value_type(1) lies in
     */
     
    /*                                                              updated: 18.12.2011 
     * a normal dice for example:
     * the constructor takes ctor(interval_start, interval_end, sections in interval)
     * 
     * NOTE: the ctor for int-types (checked with boost::is_integrale<>) is slightly
     *     different than the one for float_types (checked with boost::is_floating_point<>)
     *     This ctor gets the section length by (end - start) / (amount of sections - 1)
     *     so that the start and endpoint are valid points
     */
    
    histogram_old<int, int> a(1, 6, 6); 
    
    // put in a value via stream op
    a << 1;
    a << 2;
    a << 6;

    // put in a value via +=
    a[1] += 2;
    
    // via ++
    a[1]++;
    ++a[1];
    
    // one can also stream in a pair where the first position 
    // is the location and the second the weigth
    a << pair<int, int>(1, 1);
    
    // print
    cout << "a[1] = " << a[1] << endl;
    
    //copy ctor
    histogram_old<int, int> b(a);
    
    cout << a << endl;
    cout << b << endl;
    
    cout << "test double histo" << endl;
    cout << "-----------------" << endl;
    
    /*                                                                  updated: 18.12.2011
     * NOTE: the ctor for int-types (checked with boost::is_integrale<>) is slightly
     *     different than the one for float_types (checked with boost::is_floating_point<>)
     *     This ctor gets the section length by (end - start) / (amount of sections)
     *     the startpoint is also valid but is put in the first interval, regardless
           that a previouse interval would end there. This could mess up statistcs
     */
     
    histogram_old<double, double> fail(0, 10, 100); //calls int-ctor ==> messes up statistics
    histogram_old<double, double> c(0., 10., 100);//important to give double values to the ctor
    
    //these two values end up in the same interval
    c << 0; //startpoint
    c << 0.00000001;
    
    //same interval
    c << 1.000001;
    c << 1.099999;
    c << 1.1;
    
    
    c << 9.9001;
    c << 10; //endpoint

    // put in a value via +=
    c[1.06] += 2.1;
    
    // via ++
    c[2]++;
    ++c[2.000001];
    
    // one can also stream in a pair where the first position is 
    // the location and the second the weigth
    c << pair<double, double>(1.02, 1.4);
    
    // print
    cout << "c[1] = " << c[1] << endl;
    cout << "c[1.05] = " << c[1.000000001] << endl;
    cout << "c[1.05] = " << c[1.099999999] << endl;
    cout << "c[1.1] = " << c[1.1] << endl;
    cout << "c[0] = " << c[0] << endl;
    cout << "c[10] = " << c[10] << endl;
    
    //copy ctor
    histogram_old<double, double> d(c);
    
    cout << c << endl;
    cout << d << endl;
    
    //output:
    /*  
        test int histo
        --------------
        a[1] = 6
        histogram_old: mean: 1.75
        histogram_old: mean: 1.75
        test double histo
        -----------------
        c[1] = 5.5
        c[1.05] = 5.5
        c[1.05] = 5.5
        c[1.1] = 1
        c[0] = 2
        c[10] = 2
        histogram_old: mean: 2.5
        histogram_old: mean: 2.5
    */
}
