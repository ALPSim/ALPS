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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                             updated: 20.12.2011                 *
 * the accumulator is built modular                                                *
 * one can choose which adapters are needed                                        *
 * some adapter require others and include them automatically                      *
 *                                                                                 *
 * -----------------+-----------------------+----------------------                *
 * possible adapter | includes also         | name of the function                 *
 * -----------------+-----------------------+----------------------                *
 * Count            |                       | count(x) or x.count()                *
 * Mean             | Count                 | mean(x)  or x.mean()                 *
 * Error            | Count, Mean           | error(x) or x.error()                *
 * FixSizeBinning   | Count, Mean, Error    | fix_size_bin(x) or ..                *
 * MaxNumBinning    | Count, Mean, Error    | max_num_bin(x) or ...                *
 * LogBinning       | Count, Mean, Error    | log_bin(x) or ...                    *
 * Autocorrelation  | Count, Mean, Error    | autocorr(x)    or ...                *
 * -----------------+-----------------------+----------------------                *
 *                                                                                 *
 * NOTE: the first template parameter MUST be the value_type                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <iostream>
#include <alps/ngs.hpp>

using namespace std;
using namespace alps::alea;

int main()
{

    //is an accum with count, mean, error and fix_size_bin functions
    typedef accumulator<int, FixSizeBinning> accum; 
    
    /*                                             updated: 20.12.2011 
     * initialisation with boost::parameter
     * the tree parameter are:
     * bin_size for FixSizeBinning with default 128
     * bin_number for MaxNumBinning with default 128
     * bin_log for LogBinning with default 10
     * bin_auto for Autocorrelation with default 10
     */

    accum demo_all(bin_size = 2, bin_number = 4);
    //order with boost::parameter doesn't matter
    accum demo_all_order(bin_number = 4, bin_size = 2);
    accum demo_all_default();
    
    //ctor
    accum a(bin_size = 2);
    
    //put in values
    a << 1;
    a << 3;
    a << 5;
    a << 7;
    a << 9;
    
    //get informations
    count(a);
    mean(a);
    error(a);
    fix_size_bin(a);
    
    //print
    cout << a << endl;
    
    //copy ctor
    accum b(a);
    
    //print
    cout << b << endl;
    
    //put it in a measurement-wrapper
    measurement m(accum(bin_size = 2));
    
    //values
    m << 1;
    m << 2;
    m << 3;
    m << 4;
    
    //get informations
    m.get<int>().count();
    m.get<int>().mean();
    m.get<int>().error();
    m.get<int>().fix_size_bin();
    
    //extract accum
    m.extract<accum>();
    extract<accum>(m);
    
    //copy ctor
    measurement n(m);
    
    //print
    cout << m << endl;
    cout << n << endl;
    
    accumulator<int, Autocorrelation> v(bin_log = 1);
    
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
    v << 1;
    cout << v << endl;
}
