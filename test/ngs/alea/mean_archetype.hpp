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


#ifndef TEST_NGS_MEAN_ARCHETYPE_HEADER
#define TEST_NGS_MEAN_ARCHETYPE_HEADER

#include <alps/hdf5/archive.hpp>
#include <iostream>

using namespace std;

struct mean_archetype
{
    mean_archetype() {}
    
    mean_archetype operator+(mean_archetype const & arg)
    {
        return mean_archetype();
    }
    
    mean_archetype(int const & arg){}
    
    mean_archetype operator/(double const & arg)
    {
        return mean_archetype();
    }
    
    mean_archetype operator=(mean_archetype rhs)
    {
        return mean_archetype();
    }
    
    mean_archetype operator+=(mean_archetype rhs)
    {
        return mean_archetype();
    }

    void save(alps::hdf5::archive & ar) const {}
    void load(alps::hdf5::archive & ar) {}
};

mean_archetype operator*(mean_archetype const & arg, mean_archetype const & arg2)
{
    return mean_archetype();
}

mean_archetype operator/(mean_archetype const & arg, mean_archetype const & arg2)
{
    return mean_archetype();
}

ostream & operator<<(ostream & out, mean_archetype arg)
{
    out << "mean_archetype";
    return out;
}
#endif // TEST_NGS_MEAN_ARCHETYPE_HEADER
