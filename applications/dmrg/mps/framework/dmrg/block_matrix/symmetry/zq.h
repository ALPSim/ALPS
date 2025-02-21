/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef SYMMETRY_ZQ_H
#define SYMMETRY_ZQ_H

#include <iostream>
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


template<int Q>
class ZqCharge
{
    friend ZqCharge<Q> operator+(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return ZqCharge<Q>((a.c_+b.c_)%Q);
    }

    friend std::ostream& operator<<(std::ostream& os, ZqCharge<Q> c)
    {
        os << c.get();
        return os;
    }

    friend bool operator<(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return a.c_ < b.c_;
    }

    friend bool operator>(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return a.c_ > b.c_;
    }
    
public:
    ZqCharge(unsigned int c = 0) : c_(c) { }
    
    bool operator==(ZqCharge<Q> b) const { return c_ == b.c_; }
    bool operator!=(ZqCharge<Q> b) const { return c_ != b.c_; }
    ZqCharge operator-() const
    {
        if (c_ == 0)
            return 0;
        else
            return Q-c_;
    }
    
    unsigned int get() { return c_; }
    
    template<class Archive> void save(Archive & ar) const { ar["c"] << c_; }
    template<class Archive> void load(Archive & ar) { ar["c"] >> c_; }
    
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & c_;
    }
    
protected:
    unsigned int c_;
};

template<int Q>
class Zq
{
public:
    typedef ZqCharge<Q> charge;
    typedef int subcharge; // Used if charge is site-dependent
    
    static const charge IdentityCharge;
    static const bool finite = true;

    static const int q = Q;
    
    static charge fuse(charge a, charge b)
    {
        return a+b;
    }
    
    template<int R> static charge fuse(const boost::array<charge, R> &v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; i++)
            ret = fuse(ret, v[i]);
        return ret;
    }
};

template<int Q> const typename Zq<Q>::charge Zq<Q>::IdentityCharge = ZqCharge<Q>(0);

#endif
