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
