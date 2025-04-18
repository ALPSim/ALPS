/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MODELS_TERM_DESCRIPTOR_H
#define MODELS_TERM_DESCRIPTOR_H

#include <boost/tuple/tuple.hpp>
#include <vector>

namespace detail {
    
    struct pos_tag_lt {
        typedef boost::tuple<int, unsigned int> value_type;
        inline bool operator() (value_type const& lhs, value_type const& rhs)
        {
            return (boost::get<0>(lhs) < boost::get<0>(rhs));
        }
    };
    
}

template <typename T>
class term_descriptor : public std::vector<boost::tuple<int, unsigned int> > {
public:
    typedef int pos_type;
    typedef unsigned int tag_type;
    typedef boost::tuple<pos_type, tag_type> value_type;
    
    typedef std::vector<value_type> base;
    typedef base::size_type size_type;
    typedef typename base::iterator iterator;
    typedef typename base::const_iterator const_iterator;
    
    T coeff;
    bool is_fermionic;
    term_descriptor() : coeff(1.), is_fermionic(false) { }
    
    pos_type position(size_type i) const     { return boost::get<0>((*this)[i]); }
    tag_type operator_tag(size_type i) const { return boost::get<1>((*this)[i]); }
    
    /// utilities
    void canonical_order() // TODO: check and fix for fermions
    {
        std::sort(begin(), end(), detail::pos_tag_lt());
    }
    
    bool operator< (term_descriptor const & rhs) const
    {
        if (this->size() == 0) return true;
        if (rhs.size()   == 0) return false;
        
        if (this->position(0) == rhs.position(0))
            return this->size() > rhs.size();
        return this->position(0) < rhs.position(0);
    }
    
    bool site_match (term_descriptor const & rhs) const
    {
        if (this->size() == rhs.size())
        {
            bool ret = true;
            for (std::size_t p=0; p<this->size() && ret; ++p)
                ret = (this->position(p) == rhs.position(p));
            return ret;
        } else if (this->size() == 2 && rhs.size() == 1)
            return (this->position(0) == rhs.position(0) || this->position(1) == rhs.position(0));
        else if (this->size() == 1 && rhs.size() == 2)
            return (this->position(0) == rhs.position(0) || this->position(0) == rhs.position(1));
        else
        {
            throw std::runtime_error("site_match not implemented for this type of operator." );
            return false;
        }
        
    }
    
    bool overlap (term_descriptor const & rhs) const
    {
        return !( (boost::get<0>(*this->rbegin()) < boost::get<0>(*rhs.begin())) || (boost::get<0>(*rhs.rbegin()) < boost::get<0>(*this->begin())) );
    }
};

/// ostream
template<typename T>
std::ostream & operator<< (std::ostream & os, term_descriptor<T> const& term)
{
    os << "coeff: " << term.coeff << std::endl;
    os << "operators:";
    for (int i=0; i<term.size(); ++i)
        os << " {"  << term.position(i) << "," << term.operator_tag(i) << "}";
    os << std::endl;
    return os;
}


#endif
