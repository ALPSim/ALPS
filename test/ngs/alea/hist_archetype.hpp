/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 1997-2011 by Lukas Gamper
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
#ifndef HIST_ARCHETYPE_HEADER
#define HIST_ARCHETYPE_HEADER

#include <iostream>

using namespace std;

struct hist_archetype
{
    hist_archetype() {}
    
    operator double() const
    {
        return double();
    }
    
    hist_archetype operator+(hist_archetype const & arg) const
    {
        return hist_archetype();
    }
    hist_archetype operator-(hist_archetype const & arg) const
    {
        return hist_archetype();
    }
    
    hist_archetype(int const & arg){}
    
    hist_archetype operator/(double const & arg) const
    {
        return hist_archetype();
    }
    
    hist_archetype operator=(hist_archetype rhs)
    {
        return hist_archetype();
    }
    
    hist_archetype operator+=(hist_archetype rhs)
    {
        return hist_archetype();
    }
};
hist_archetype operator*(hist_archetype const & arg, hist_archetype const & arg2)
{
    return hist_archetype();
}
hist_archetype operator*(hist_archetype const & arg, unsigned int const & arg2)
{
    return hist_archetype();
}

hist_archetype operator*(int const & arg, hist_archetype const & arg2)
{
    return hist_archetype();
}

ostream & operator<<(ostream & out, hist_archetype arg)
{
    out << "hist_archetype";
    return out;
}
#endif // HIST_ARCHETYPE_HEADER
