/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MEASUREMENT_TERM_DESC_HPP
#define MEASUREMENT_TERM_DESC_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

template <class T>
struct MeasurementTermDescriptor {
    typedef std::vector<std::string> op_prod_type;
    T coeff;
    std::vector<op_prod_type> op_names;
    
    MeasurementTermDescriptor() { }
    MeasurementTermDescriptor(T const& coeff_, std::vector<op_prod_type> const& op_names_)
    : coeff(coeff_)
    , op_names(op_names_)
    { }
};

template <class T>
std::ostream & operator<<(std::ostream & os, MeasurementTermDescriptor<T> const& meas_term_desc)
{
    os << "coeff=" << meas_term_desc.coeff;
    os << ", op_names=[";
    for (unsigned i=0; i<meas_term_desc.op_names.size(); ++i) {
        os << "[";
        std::copy(meas_term_desc.op_names[i].begin(), meas_term_desc.op_names[i].end()-1, std::ostream_iterator<std::string const&>(os, ","));
        os << meas_term_desc.op_names[i].back();
        os << "]";
        if (i < meas_term_desc.op_names.size()-1)
            os << ", ";
    }
    os << "]";
    return os;
}

template <class T>
std::ostream & operator<<(std::ostream & os, std::vector<MeasurementTermDescriptor<T> > const& meas_terms_desc)
{
    os << "[\n  ";
    std::copy(meas_terms_desc.begin(), meas_terms_desc.end()-1, std::ostream_iterator<MeasurementTermDescriptor<T> const&>(os, ",\n  "));
    os << meas_terms_desc.back();
    os << "\n]\n";
    return os;
}

#endif