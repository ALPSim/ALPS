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
 *                                             updated: 21.05.2012                 *
 * the accumulator is built modular                                                *
 * one can choose which implementations are needed                                 *
 * some implementation require others and include them automatically               *
 *                                                                                 *
 * -----------------+-----------------------+----------------------                *
 * possible impl    | includes also         | name of the function                 *
 * -----------------+-----------------------+----------------------                *
 * mean             | Count                 | mean(x)                              *
 * error            | Count, mean           | error(x)                             *
 * fixed_size_binn  | Count, mean, error    | fixed_size_binning(x)                *
 * MaxNumBinning    | Count, mean, error    | max_num_binning(x)                   *
 * log_binning      | Count, mean, error    | log_binning(x)                       *
 * autocorrelation  | Count, mean, error    | autocorrelation(x)                   *
 * -----------------+-----------------------+----------------------                *
 *                                                                                 *
 * NOTE: the first template parameter MUST be the value_type                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <iostream>
#include <alps/ngs.hpp>

using namespace std;
using namespace alps::accumulator;

int main()
{
    
    //accumulator with all the features so far
    typedef accumulator<  int
                        , features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                //~ , tag::converged
                                //~ , tag::tau
                                , tag::histogram
                                >
                        > accum;
    
    /*                                             updated: 21.05.12 
     * initialisation with boost::parameter
     * the two parameter are:
     * bin_size for tag::fixed_size_binning with default 128
     * bin_num for MaxNumBinning with default 128
     */
    
    //------------------- constructor -------------------
    accum demo_all_default;
    accum demo_partial_1(bin_size = 2);     //bin_num is default 128
    accum demo_partial_2(bin_num = 2);   //bin_size is default 128
    accum demo(bin_size = 2, bin_num = 2);
    
    //order with boost::parameter doesn't matter
    accum demo_all_order(bin_num = 4, bin_size = 2);
    
    //------------------- stream operator -------------------
    demo << 1;
    demo << 3;
    demo << 5;
    demo << 7;
    demo << 9;
    
    //------------------- get informations -------------------
    cout << count(demo) << endl; 
    //output:   5
    cout << mean(demo) << endl;
    //output:   5
    cout << error(demo) << endl;
    //output:   1.41421
    cout << fixed_size_binning(demo).bins()[0] << endl;
    //output:   2     //(1+3)/2
    cout << max_num_binning(demo).bins()[0] << endl;
    //output:   2     //(1+3)/2
    cout << log_binning(demo).bins()[0] << endl;
    //output:   1     //first element
    cout << autocorrelation(demo).error(2) << endl;
    cout << autocorrelation(demo).bins()[1] << endl;
    cout << autocorrelation(demo).sum()[1] << endl;
    
    cout << converged(demo) << endl;
    //output:   2 (alps::accumulator::error_convergence maybe)   //not jet implemented
    cout << tau(demo) << endl;
    //output:   42    //not jet implemented
    cout << histogram(demo) << endl;
    //output:   272.15    //not jet implemented
    
    //------------------- copy constructor -------------------
    accum copy(demo);
    
    //------------------- print to an ostream -------------------
    cout << demo << endl;
    //output:   ValueType: i Count: 5 tag::mean: 5 tag::error: 1.41421 FixBinSize->BinNumber: 2 
    //          MaxBinningNumber->MaxBinNumber: 2 Log Binning: tag::autocorrelation:
    cout << copy << endl;
    //output:   ValueType: i Count: 5 tag::mean: 5 tag::error: 1.41421 FixBinSize->BinNumber: 2 
    //          MaxBinningNumber->MaxBinNumber: 2 Log Binning: tag::autocorrelation:
    
    
    //------------------- detail::accumulator_wrapper -------------------
    //put an accumulator in a detail::accumulator_wrapper-wrapper (copy)
    detail::accumulator_wrapper m_from_accum(demo);
    //construct in initialization
    detail::accumulator_wrapper m(accum(bin_size = 4, bin_num = 10));
    
    //------------------- stream operator -------------------
    for(int i = 0; i < 64; ++i)
        m << 1;
    
    //------------------- get informations -------------------
    cout << count(m.get<int>()) << endl;
    //output:   64
    cout << mean(m.get<int>()) << endl;
    //output:   1
    cout << error(m.get<int>()) << endl;
    //output:   0
    cout << tau(m.get<int>()) << endl;
    cout << histogram(m.get<int>()) << endl;
    
    //------------------- extract accum -------------------
    accum extraction(m.extract<accum>());
    accum extraction_free(extract<accum>(m));
    
    //------------------- copy constructor -------------------
    detail::accumulator_wrapper n(m);
    
    
    //------------------- print -------------------
    cout << m << endl;
    //output:   ValueType: i Count: 64 tag::mean: 1 tag::error: 0 FixBinSize->BinNumber: 16 
    //          MaxBinningNumber->MaxBinNumber: 10 Log Binning: tag::autocorrelation: 
    cout << n << endl;
    //output:   ValueType: i Count: 64 tag::mean: 1 tag::error: 0 FixBinSize->BinNumber: 16 
    //          MaxBinningNumber->MaxBinNumber: 10 Log Binning: tag::autocorrelation: 
    
    //------------------- binnings -------------------    
    
    cout << n.count() << endl;
    cout << fixed_size_binning(n.get<int>()).bins()[0] << " ";
    cout << n.get<int>().fixed_size_binning().bins()[1] << " ";
    cout << n.get<int>().fixed_size_binning().bins()[2] << " ";
    cout << n.get<int>().fixed_size_binning().bin_size() << endl;
    //~ //output: 1 1 1 1
    cout << n.get<int>().max_num_binning().bins()[0] << " ";
    cout << n.get<int>().max_num_binning().bins()[1] << " ";
    cout << n.get<int>().max_num_binning().bins()[2] << " ";
    cout << n.get<int>().max_num_binning().bins()[3] << endl;
    //~ //output: 1 1 1 1
    cout << n.get<int>().log_binning().bins()[0] << " ";
    cout << n.get<int>().log_binning().bins()[1] << " ";
    cout << n.get<int>().log_binning().bins()[2] << " ";
    cout << n.get<int>().log_binning().bins()[3] << endl;
    //~ //output: 1 1 1 1
    cout << n.get<int>().autocorrelation() << " ";
    cout << n.get<int>().autocorrelation() << " ";
    cout << n.get<int>().autocorrelation() << " ";
    cout << n.get<int>().autocorrelation() << endl;
    //output: 64 128 256 512

    //------------------- doesn't work jet -------------------
    //~ //need to adjust resulttype-wrapper
    //~ cout << converged(n.get<int>()) << endl;
    //~ cout << tau(n.get<int>()) << endl;
    
    
    int i = 0;
    
    cout << i << endl;
    
    //~ detail::weak_type_ptr a(i);
    //~ 
    //~ detail::make_data<int> b(5);
    //~ 
    //~ ++b().cast<int>();
    //~ 
    //~ cout << b.get() << endl;
    //~ 
    //~ ++a.cast<int>();
    //~ 
    //~ cout << i << endl;
    
    accumulator<int> aa;
}
