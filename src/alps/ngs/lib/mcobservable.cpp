/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/mcobservable.hpp>

#include <alps/alea/observable.h>

#include <vector>
#include <valarray>
#include <iostream>

namespace alps {

    mcobservable::mcobservable()
        : impl_(NULL) 
    {}

    mcobservable::mcobservable(Observable const * obs) {
        ref_cnt_[impl_ = obs->clone()] = 1;
    }

    mcobservable::mcobservable(mcobservable const & rhs) {
        ++ref_cnt_[impl_ = rhs.impl_];
    }

    mcobservable::~mcobservable() {
        if (impl_ && !--ref_cnt_[impl_])
            delete impl_;
    }

    mcobservable & mcobservable::operator=(mcobservable rhs) {
        if (impl_ && !--ref_cnt_[impl_])
            delete impl_;
        ++ref_cnt_[impl_ = rhs.impl_];
        return *this;
    }

    Observable * mcobservable::get_impl() {
        return impl_;
    }

    Observable const * mcobservable::get_impl() const {
        return impl_;
    }

    std::string const & mcobservable::name() const {
        return impl_->name();
    }

     template<> ALPS_DECL mcobservable & mcobservable::operator<< <double>(double const & value) {
        (*impl_) << value;
        return *this;
    }

     template<> ALPS_DECL mcobservable & mcobservable::operator<< <std::vector<double> >(std::vector<double>  const & value) {
        std::valarray<double> varr(value.size());
        std::copy(value.begin(), value.end(), &varr[0]);
        (*impl_) << varr;
        return *this;
    }

     template<> ALPS_DECL mcobservable & mcobservable::operator<< <std::valarray<double> >(std::valarray<double>  const & value) {
        (*impl_) << value;
        return *this;
    }

    void mcobservable::save(hdf5::archive & ar) const {
        impl_->save(ar);
    }

    void mcobservable::load(hdf5::archive & ar) {
        impl_->load(ar);
    }

    void mcobservable::merge(mcobservable const & obs) {
        if (!impl_->can_merge()) {
            Observable* unmergeable = impl_;
            ++ref_cnt_[impl_ = unmergeable->convert_mergeable()];
            if (!--ref_cnt_[unmergeable])
                delete unmergeable;
        }
        impl_->merge(*obs.get_impl());
    }

    void mcobservable::output(std::ostream & os) const {
        os << *(impl_);
    }

    std::map<Observable *, std::size_t> mcobservable::ref_cnt_;

    std::ostream & operator<<(std::ostream & os, mcobservable const & obs) {
        obs.output(os);
        return os;
    }

}
