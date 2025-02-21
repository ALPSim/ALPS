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


#ifndef TEST_NGS_ERROR_ARCHETYPE_HEADER
#define TEST_NGS_ERROR_ARCHETYPE_HEADER

#include "mean_archetype.hpp"

#include <iostream>

using namespace std;

struct error_archetype
{
    error_archetype() {}
    
    error_archetype operator+(error_archetype const & arg)
    {
        return error_archetype();
    }
    
    error_archetype(int const & arg){}
    
    error_archetype operator/(double const & arg) const
    {
        return error_archetype();
    }
    
    error_archetype operator=(error_archetype rhs)
    {
        return error_archetype();
    }
    
    error_archetype operator+=(error_archetype rhs)
    {
        return error_archetype();
    }

    void save(alps::hdf5::archive & ar) const {}
    void load(alps::hdf5::archive & ar) {}
};

error_archetype operator*(error_archetype const & arg, error_archetype const & arg2)
{
    return error_archetype();
}

error_archetype operator/(error_archetype const & arg, error_archetype const & arg2)
{
    return error_archetype();
}

error_archetype operator-(error_archetype const & arg, error_archetype const & arg2)
    {
        return error_archetype();
    }

ostream & operator<<(ostream & out, error_archetype arg)
{
    out << "error_archetype";
    return out;
}

error_archetype sqrt(error_archetype)
{
    return error_archetype();
}
#endif // TEST_NGS_ERROR_ARCHETYPE_HEADER
