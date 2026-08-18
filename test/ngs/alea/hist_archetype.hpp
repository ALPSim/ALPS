/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 1997-2011 by Lukas Gamper
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
