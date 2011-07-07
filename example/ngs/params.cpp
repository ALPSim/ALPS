/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/ngs/hdf5.hpp>
#include <alps/ngs/params.hpp>

#include <fstream>
#include <iostream>
#include <iterator>

int main(int argc, char *argv[]) {

    // load a hdf5 file into a parameter object
    alps::hdf5::archive ar("param.h5");
    alps::params params_h5(ar);

    std::vector<std::string> keys_h5 = params_h5.keys();
    for (std::vector<std::string>::const_iterator it = keys_h5.begin(); it != keys_h5.end(); ++it)
        std::cout << *it << ":" << params_h5[*it] << std::endl;

    std::cout << std::endl;

    // load a text file into a parameter object
    std::ifstream file;
    file.open("param.txt");
    std::string txt;
    while(file.good()) {
        std::string line;
        file >> line;
        txt += line + "\n";
    }
    file.close();
    alps::params params_txt(txt);

    std::vector<std::string> keys_txt = params_txt.keys();
    for (std::vector<std::string>::const_iterator it = keys_txt.begin(); it != keys_txt.end(); ++it)
        std::cout << *it << ":" << params_txt[*it] << std::endl;

}