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

#ifndef SYMMETRY_Z2_H
#define SYMMETRY_Z2_H

#include <iostream>

#include <boost/functional/hash.hpp>
 
#include <alps/hdf5.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

class Ztwo {
	public:
		typedef enum { Plus = 0, Minus = 1 } charge;
        typedef int subcharge; // used if charge is site-dependent
		
		static const charge IdentityCharge = Plus;
        static const bool finite = true;
		
		static inline charge fuse(charge a, charge b)
		{
			if (a == b)
				return Plus;
			else
				return Minus;
		}
		
		template<int R>
		static charge fuse(boost::array<charge, R> v)
		{
			// this operation actually could be rearranged into a tree
			for (int i = 1; i < R; i++)
				v[0] = fuse(v[0], v[i]);
			return v[0];
		}
};

inline void save(alps::hdf5::archive & ar,
                 std::string const & p,
                 Ztwo::charge const & v,
                 std::vector<std::size_t> size = std::vector<std::size_t>(),
                 std::vector<std::size_t> chunk = std::vector<std::size_t>(),
                 std::vector<std::size_t> offset = std::vector<std::size_t>())
{
    ar[p] << static_cast<int>(v);
}

inline void load(alps::hdf5::archive & ar,
                 std::string const & p,
                 Ztwo::charge & v,
                 std::vector<std::size_t> size = std::vector<std::size_t>(),
                 std::vector<std::size_t> chunk = std::vector<std::size_t>(),
                 std::vector<std::size_t> offset = std::vector<std::size_t>())
{
    int t;
    ar[p] >> t;
    v = (t == 0 ? Ztwo::Plus : Ztwo::Minus);
}

template <class Archive>
inline void serialize(Archive & ar, Ztwo::charge & c, const unsigned int version)
{
    ar & c;
}

inline Ztwo::charge operator-(Ztwo::charge a) { return a; }

inline std::ostream& operator<<(std::ostream& ost, Ztwo::charge c)
{
	if (c == Ztwo::Plus)
		ost << "Plus";
	else if (c == Ztwo::Minus)
		ost << "Minus";
	else
		ost << "???";
	return ost;
}
inline std::ostream& operator<<(std::ostream& ost, const std::vector<Ztwo::charge> &c)
{
	ost << "[ ";
	for (std::vector<Ztwo::charge>::const_iterator it = c.begin();
		 it != c.end();
		 it++)
		ost << ", " << *it;
	ost << " ]";
	return ost;
}

#endif
