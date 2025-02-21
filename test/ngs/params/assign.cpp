/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
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

#include <alps/ngs/params.hpp>

int main() {

    alps::params parms;
    parms["char"] = static_cast<char>(1);
    parms["signed char"] = static_cast<signed char>(1);
    parms["unsigned char"] = static_cast<unsigned char>(1);
    parms["short"] = static_cast<short>(1);
    parms["unsigned short"] = static_cast<unsigned short>(1);
    parms["int"] = static_cast<int>(1);
    parms["unsigned"] = static_cast<unsigned>(1);
    parms["long"] = static_cast<long>(1);
    parms["unsigned long"] = static_cast<unsigned long>(1);
    parms["long long"] = static_cast<long long>(1);
    parms["unsigned long long"] = static_cast<unsigned long long>(1);
    parms["float"] = static_cast<float>(1);
    parms["double"] = static_cast<double>(1);
    parms["long double"] = static_cast<long double>(1);
    parms["bool"] = static_cast<bool>(1);
    parms["std::string"] = std::string("asdf");

    std::cout << parms << std::endl;
    return 0;
}
