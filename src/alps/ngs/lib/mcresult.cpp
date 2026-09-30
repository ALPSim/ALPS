/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/mcresult.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/ngs/numeric/array.hpp>
#include <alps/ngs/lib/mcresult_impl_base.ipp>
#include <alps/ngs/lib/mcresult_impl_derived.ipp>

// #ifdef ALPS_NGS_USE_NEW_ALEA
//     #include <alps/ngs/alea/accumulator_set.hpp>
// #endif

#include <alps/alea/observable.h>
#include <alps/alea/abstractsimpleobservable.h>

#include <iostream>
#include <stdexcept>
#include <memory>
#include <utility>

namespace alps {

    detail::mcresult_impl_base * mcresult::implementation() const {
        if (!impl_)
            throw std::runtime_error("Result has no measurements" + ALPS_STACKTRACE);
        return impl_;
    }

    mcresult::mcresult()
        : impl_(NULL) 
    {}

    mcresult::mcresult(Observable const * obs) {
        construct(obs);
    }
    
// #ifdef ALPS_NGS_USE_NEW_ALEA
//     template<typename T> bool has_cast(alps::accumulator::detail::accumulator_wrapper const & acc_wrapper) { 
//         try { 
//             acc_wrapper.get<T>();
//             return true;
//         }  catch(std::bad_cast)  { 
//             return false;
//         }
//     }
    
//     mcresult::mcresult(alps::accumulator::detail::accumulator_wrapper const & acc_wrapper) {
//         //------------------- find out value_type  -------------------
//         if(has_cast<double>(acc_wrapper))
//             impl_ = new detail::mcresult_impl_derived<detail::mcresult_impl_base, double>(acc_wrapper);
//         else if(has_cast<std::vector<double> >(acc_wrapper))
//            impl_ = new detail::mcresult_impl_derived<detail::mcresult_impl_base, std::vector<double> >(acc_wrapper);
// /*
//         else if(has_cast<boost::array<double, 3> >(acc_wrapper))
//             impl_ = new detail::mcresult_impl_derived<detail::mcresult_impl_base, boost::array<double, 3> >(acc_wrapper);
// */
//         else
//             throw std::runtime_error("type not supported in mcresult-constructor" + ALPS_STACKTRACE);
//         ref_cnt_[impl_] = 1;
//     }
// #endif

    mcresult::mcresult(mcresult const & rhs) {
        impl_ = rhs.impl_;
        if (impl_) ++ref_cnt_[impl_];
    }

    mcresult::mcresult(mcobservable const & obs) {
        construct(obs.get_impl());
    }

    mcresult::~mcresult() {
        if (impl_ && !--ref_cnt_[impl_]) {
            ref_cnt_.erase(impl_);
            delete impl_;
        }
    }

    mcresult & mcresult::operator=(mcresult rhs) {
        std::swap(impl_, rhs.impl_);
        return *this;
    }
    
    #define ALPS_MCRESULT_TPL_IMPL(T)                                                                                               \
        template<> ALPS_DECL bool mcresult::is_type< T >() const { return impl_ && implementation()->is_type< T >(); }                         \
        template<> ALPS_DECL std::vector< T > const & mcresult::bins< T >() const { return implementation()->bins< T >(); }                    \
        template<> ALPS_DECL T const & mcresult::mean< T >() const { return implementation()->mean< T >(); }                                   \
        template<> ALPS_DECL T const & mcresult::error< T >() const { return implementation()->error< T >(); }                                 \
        template<> ALPS_DECL T const & mcresult::variance< T >() const { return implementation()->variance< T >(); }                           \
        template<> ALPS_DECL T const & mcresult::tau< T >() const { return implementation()->tau< T >(); }                                     \
        template<> ALPS_DECL covariance_type<T>::type mcresult::covariance< T >(mcresult const & arg) const {                       \
            return implementation()->covariance< T >(*arg.implementation());                                                                              \
        }                                                        \
        template<> ALPS_DECL covariance_type<T>::type mcresult::accurate_covariance< T >(mcresult const & arg) const {                       \
            return implementation()->accurate_covariance< T >(*arg.implementation());                                                                           \
        }                                                        \
        template<> ALPS_DECL mcresult & mcresult::operator+=< T >( T const & rhs) { implementation()->add_assign(rhs); return *this; }         \
        template<> ALPS_DECL mcresult & mcresult::operator-=< T >( T const & rhs) { implementation()->sub_assign(rhs); return *this; }         \
        template<> ALPS_DECL mcresult & mcresult::operator*=< T >( T const & rhs) { implementation()->mul_assign(rhs); return *this; }         \
        template<> ALPS_DECL mcresult & mcresult::operator/=< T >( T const & rhs) { implementation()->div_assign(rhs); return *this; }
    ALPS_MCRESULT_TPL_IMPL(double)
    ALPS_MCRESULT_TPL_IMPL(std::vector<double>)
    #undef ALPS_MCRESULT_TPL_IMPL

    #define ALPS_NGS_MCRESULT_OPERATOR_IMPL(OP, NAME)         \
        mcresult & mcresult:: OP (mcresult const & rhs) {     \
            implementation()-> NAME ## _assign (rhs.implementation());              \
            return *this;                                     \
        }
    ALPS_NGS_MCRESULT_OPERATOR_IMPL(operator+=, add)
    ALPS_NGS_MCRESULT_OPERATOR_IMPL(operator-=, sub)
    ALPS_NGS_MCRESULT_OPERATOR_IMPL(operator*=, mul)
    ALPS_NGS_MCRESULT_OPERATOR_IMPL(operator/=, div)
    #undef ALPS_NGS_MCRESULT_OPERATOR_IMPL

    bool mcresult::can_rebin() const {
        return implementation()->can_rebin();
    }

    bool mcresult::jackknife_valid() const {
        return implementation()->jackknife_valid();
    }

    uint64_t mcresult::count() const {
        return impl_ ? implementation()->count() : 0;
    }

    uint64_t mcresult::bin_size() const {
        return implementation()->bin_size();
    }

    uint64_t mcresult::max_bin_number() const {
        return implementation()->max_bin_number();
    }

    std::size_t mcresult::bin_number() const {
        return implementation()->bin_number();
    }

    bool mcresult::has_variance() const {
        return implementation()->has_variance();
    }

    bool mcresult::has_tau() const {
        return implementation()->has_tau();
    }

    void mcresult::set_bin_size(uint64_t binsize) {
        implementation()->set_bin_size(binsize);
    }

    void mcresult::set_bin_number(uint64_t bin_number) {
        implementation()->set_bin_number(bin_number);
    }

    void mcresult::save(hdf5::archive & ar) const {
        if (!impl_)
            throw std::runtime_error("Cannot save an uninitialized result" + ALPS_STACKTRACE);
        implementation()->save(ar);
    }

    void mcresult::load(hdf5::archive & ar) {
        // Read into an independent payload before touching this result. This
        // also discovers the type of a default-constructed result and keeps
        // aliases and the old value intact if any archive read fails.
        std::unique_ptr<detail::mcresult_impl_base> payload;
        if (ar.is_scalar("mean/value"))
            payload.reset(new detail::mcresult_impl_derived<detail::mcresult_impl_base, double>(alea::mcdata<double>()));
        else if (ar.dimensions("mean/value") == 1)
            payload.reset(new detail::mcresult_impl_derived<detail::mcresult_impl_base, std::vector<double> >(alea::mcdata<std::vector<double> >()));
        else
            throw std::runtime_error("Unsupported result shape" + ALPS_STACKTRACE);
        payload->load(ar);
        mcresult replacement;
        ref_cnt_[payload.get()] = 1;
        replacement.impl_ = payload.release();
        std::swap(impl_, replacement.impl_);
    }

    void mcresult::output(std::ostream & os) const {
        implementation()->output(os);
    }

    #ifdef ALPS_HAVE_MPI
        mcresult mcresult::reduce(boost::mpi::communicator const & communicator, std::size_t binnumber) {
            mcresult lhs;
            detail::mcresult_impl_base * impl = implementation()->reduce(communicator, binnumber);
            if (communicator.rank() == 0)
                ref_cnt_[lhs.impl_ = impl] = 1;
            return lhs;
        }
    #endif

    bool mcresult::operator== (mcresult const & rhs) const {
        return implementation()->operator== (rhs.implementation());
    }
    bool mcresult::operator!= (mcresult const & rhs) const {
        return implementation()->operator!= (rhs.implementation());
    }

    mcresult & mcresult::operator+() {
        return *this;
    }
    mcresult & mcresult::operator-() {
        implementation()->operator-();
        return *this;
    }

    void mcresult::construct(Observable const * obs) {
        if (dynamic_cast<AbstractSimpleObservable<double> const *>(obs) != NULL)
            impl_ = new detail::mcresult_impl_derived<detail::mcresult_impl_base, double>(
                dynamic_cast<AbstractSimpleObservable<double> const &>(*obs)
            );
        else if (dynamic_cast<AbstractSimpleObservable<std::valarray<double> > const *>(obs) != NULL)
            impl_ = new detail::mcresult_impl_derived<detail::mcresult_impl_base, std::vector<double> >(
                dynamic_cast<AbstractSimpleObservable<std::valarray<double> > const &>(*obs)
            );
        else
            throw std::runtime_error("unknown observable type" + ALPS_STACKTRACE);
        ref_cnt_[impl_] = 1;
    }
    

    std::map<detail::mcresult_impl_base *, std::size_t> mcresult::ref_cnt_;

    std::ostream & operator<<(std::ostream & os, mcresult const & res) {
        res.output(os);
        return os;
    }

    #define ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(FUN_NAME)                 \
        mcresult FUN_NAME (mcresult rhs) {                               \
            mcresult lhs;                                                \
            lhs.ref_cnt_[lhs.impl_ = rhs.implementation()-> FUN_NAME ()] = 1;       \
            return lhs;                                                  \
        }
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(sin)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(cos)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(tan)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(sinh)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(cosh)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(tanh)
    // ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(asin)
    // ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(acos)
    // ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(atan)
// asinh, aconsh and atanh are not part of C++03 standard
//    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(asinh)
//    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(acosh)
//    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(atanh)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(abs)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(sq)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(cb)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(sqrt)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(cbrt)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(exp)
    ALPS_NGS_MCRESULT_FREE_UNITARY_FUN(log)
    #undef ALPS_NGS_MCRESULT_FREE_UNITARY_FUN

    mcresult pow(mcresult rhs, double exponent) {
        mcresult lhs;
        lhs.ref_cnt_[lhs.impl_ = rhs.implementation()->pow(exponent)] = 1;
        return lhs;
    }

    #define ALPS_NGS_MCRESULT_FREE_OPERATOR_TPL_IMPL(T, OP, NAME)                  \
        mcresult OP(mcresult const & lhs, T const & rhs) {                         \
            mcresult res;                                                          \
            res.ref_cnt_[res.impl_ = lhs.implementation()-> NAME (rhs)] = 1;                  \
            return res;                                                            \
        }                                                                          \
        mcresult OP(T const & lhs, mcresult const & rhs) {                         \
            mcresult res;                                                          \
            res.ref_cnt_[res.impl_ = rhs.implementation()-> NAME ## _inverse (lhs)] = 1;      \
            return res;                                                            \
        }
    #define ALPS_NGS_MCRESULT_FREE_OPERATOR_IMPL(OP, NAME)                         \
        ALPS_NGS_MCRESULT_FREE_OPERATOR_TPL_IMPL(double, OP, NAME)                 \
        ALPS_NGS_MCRESULT_FREE_OPERATOR_TPL_IMPL(std::vector<double>, OP, NAME)    \
        mcresult OP (mcresult const & lhs, mcresult const & rhs) {                 \
            mcresult res;                                                          \
            res.ref_cnt_[res.impl_ = lhs.implementation()-> NAME (rhs.implementation())] = 1;            \
            return res;                                                            \
        }
    ALPS_NGS_MCRESULT_FREE_OPERATOR_IMPL(operator+, add)
    ALPS_NGS_MCRESULT_FREE_OPERATOR_IMPL(operator-, sub)
    ALPS_NGS_MCRESULT_FREE_OPERATOR_IMPL(operator*, mul)
    ALPS_NGS_MCRESULT_FREE_OPERATOR_IMPL(operator/, div)
    #undef ALPS_MCRESULT_OPERATOR_IMPL
    #undef ALPS_NGS_MCRESULT_FREE_OPERATOR_TPL_IMPL

}
