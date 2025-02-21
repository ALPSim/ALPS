/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef LATTICE_H
#define LATTICE_H

#include "dmrg/utils/BaseParameters.h"

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

/// lattice common base
class lattice_impl {
public:
    typedef int pos_t;
    
    virtual ~lattice_impl() {}
    
    virtual std::vector<pos_t> forward(pos_t) const = 0;
    virtual std::vector<pos_t> all(pos_t) const = 0;
    
    // non-virtual!
    template<class T> T get_prop(std::string property,
                                 pos_t site) const
    {
        return boost::any_cast<T>(get_prop_(property, std::vector<pos_t>(1, site)));
    }
    
    template<class T> T get_prop(std::string property,
                                 pos_t bond1, pos_t bond2) const
    {
        std::vector<pos_t> v(2);
        v[0] = bond1; v[1] = bond2;
        return boost::any_cast<T>(get_prop_(property, v));
    }
    
    template<class T> T get_prop(std::string property,
                                 std::vector<pos_t> const & positions) const
    {
        return boost::any_cast<T>(get_prop_(property, positions));
    }
    
    // virtual!
    virtual boost::any get_prop_(std::string const &, std::vector<pos_t> const &) const = 0;
    
    virtual pos_t size() const = 0;
    virtual int maximum_vertex_type() const = 0;

};


/// lattice factory
boost::shared_ptr<lattice_impl>
lattice_factory(BaseParameters & parms);


/// pimpl resolved Lattice
class Lattice {
    typedef lattice_impl impl_type;
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
public:
    typedef impl_type::pos_t pos_t;
    
    Lattice() { }
    
    Lattice(BaseParameters & parms)
    : impl_(lattice_factory(parms))
    { }
    
    Lattice(impl_ptr impl) : impl_(impl) { }
    
    impl_ptr impl() const { return impl_; }
    
    std::vector<pos_t> forward(pos_t site) const { return impl_->forward(site); }
    std::vector<pos_t> all(pos_t site) const { return impl_->all(site); }
    
    template<class T> T get_prop(std::string property, pos_t site) const
    { return impl_->get_prop<T>(property, site); }
    
    template<class T> T get_prop(std::string property, pos_t bond1, pos_t bond2) const
    { return impl_->get_prop<T>(property, bond1, bond2); }
    
    template<class T> T get_prop(std::string property, std::vector<pos_t> const & positions) const
    { return impl_->get_prop<T>(property, positions); }
    
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & positions) const
    { return impl_->get_prop_(property, positions); }
    
    pos_t size() const { return impl_->size(); }
    int maximum_vertex_type() const { return impl_->maximum_vertex_type(); };

    
private:
    impl_ptr impl_;
};

#endif
