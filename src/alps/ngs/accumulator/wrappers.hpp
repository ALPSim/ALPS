/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_ACCUMULATOR_WRAPPER_HPP
#define ALPS_NGS_ACCUMULATOR_WRAPPER_HPP

#include <alps/ngs/accumulator/feature/mean.hpp>
#include <alps/ngs/accumulator/feature/error.hpp>
#include <alps/ngs/accumulator/feature/count.hpp>
#include <alps/ngs/accumulator/feature/weight.hpp>
#include <alps/ngs/accumulator/feature/max_num_binning.hpp>
#include <alps/ngs/accumulator/feature/binning_analysis.hpp>

#include <alps/hdf5/archive.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/variant/variant.hpp>

#include <typeinfo>
#include <stdexcept>

namespace alps {
    namespace accumulator {

        template<typename A> class derived_wrapper;

        namespace detail {
            template<typename T> struct value_wrapper {
                typedef T value_type;
            };

            typedef boost::make_variant_over<
                boost::mpl::vector<ALPS_ACCUMULATOR_VALUE_TYPES>
            >::type weight_variant_type;
        }

        template<typename T> class base_wrapper : public 
            impl::BaseWrapper<T, weight_tag, 
            impl::BaseWrapper<T, max_num_binning_tag, 
            impl::BaseWrapper<T, binning_analysis_tag, 
            impl::BaseWrapper<T, error_tag, 
            impl::BaseWrapper<T, mean_tag, 
            impl::BaseWrapper<T, count_tag, 
            detail::value_wrapper<T>
        > > > > > > {

            public:
                using typename detail::value_wrapper<T>::value_type;

                virtual ~base_wrapper() {}

                virtual void operator()(value_type const & value) = 0;
                virtual void operator()(value_type const & value, detail::weight_variant_type const & weight) = 0;

                virtual void save(hdf5::archive & ar) const = 0;
                virtual void load(hdf5::archive & ar) = 0;

                virtual void print(std::ostream & os) const = 0;
                virtual void reset() = 0;

#ifdef ALPS_HAVE_MPI
                virtual void collective_merge(boost::mpi::communicator const & comm, int root) = 0;
#endif

                virtual base_wrapper * clone() const = 0;
                virtual base_wrapper * result() const = 0;

                template<typename A> A & extract() {
                    return dynamic_cast<derived_wrapper<A> &>(*this).extract();
                }
                template<typename A> A const & extract() const {
                    return dynamic_cast<derived_wrapper<A> const &>(*this).extract();
                }

                virtual void operator+=(base_wrapper const &) = 0;
                virtual void operator-=(base_wrapper const &) = 0;
                virtual void operator*=(base_wrapper const &) = 0;
                virtual void operator/=(base_wrapper const &) = 0;

                virtual void operator+=(double) = 0;
                virtual void operator-=(double) = 0;
                virtual void operator*=(double) = 0;
                virtual void operator/=(double) = 0;

                virtual void negate() = 0;
                virtual void inverse() = 0;

                virtual void sin() = 0;
                virtual void cos() = 0;
                virtual void tan() = 0;
                virtual void sinh() = 0;
                virtual void cosh() = 0;
                virtual void tanh() = 0;
                virtual void asin() = 0;
                virtual void acos() = 0;
                virtual void atan() = 0;
                virtual void abs() = 0;
                virtual void sqrt() = 0;
                virtual void log() = 0;
                virtual void sq() = 0;
                virtual void cb() = 0;
                virtual void cbrt() = 0;
        };

        namespace detail {
            template<typename A> class foundation_wrapper : public base_wrapper<typename value_type<A>::type> {

                public:
                    foundation_wrapper(A const & arg): m_data(arg) {}

                protected:
                    A m_data;
            };
        }

        template<typename A> class derived_wrapper : public 
            impl::DerivedWrapper<A, weight_tag, 
            impl::DerivedWrapper<A, max_num_binning_tag, 
            impl::DerivedWrapper<A, binning_analysis_tag, 
            impl::DerivedWrapper<A, error_tag, 
            impl::DerivedWrapper<A, mean_tag, 
            impl::DerivedWrapper<A, count_tag, 
        detail::foundation_wrapper<A> > > > > > > {

            using typename detail::value_wrapper<typename value_type<A>::type>::value_type;
            
            public:
                derived_wrapper()
                    : 
                        impl::DerivedWrapper<A, weight_tag, 
                        impl::DerivedWrapper<A, max_num_binning_tag, 
                        impl::DerivedWrapper<A, binning_analysis_tag, 
                        impl::DerivedWrapper<A, error_tag, 
                        impl::DerivedWrapper<A, mean_tag, 
                        impl::DerivedWrapper<A, count_tag, 
                    detail::foundation_wrapper<A> > > > > > >() 
                {}

                derived_wrapper(A const & arg)
                    : 
                        impl::DerivedWrapper<A, weight_tag, 
                        impl::DerivedWrapper<A, max_num_binning_tag, 
                        impl::DerivedWrapper<A, binning_analysis_tag, 
                        impl::DerivedWrapper<A, error_tag, 
                        impl::DerivedWrapper<A, mean_tag, 
                        impl::DerivedWrapper<A, count_tag, 
                    detail::foundation_wrapper<A> > > > > > >(arg) 
                {}

                A & extract() {
                    return this->m_data;
                }
                A const & extract() const {
                    return this->m_data;
                }

                void operator()(value_type const & value) {
                    this->m_data(value);
                }

            private:
                template<typename Q>
                struct call_2_visitor: public boost::static_visitor<> {
                    call_2_visitor(A & d, value_type const & v) : data(d), value(v) {}
                    template<typename X> void operator()(X const & arg) const {
                        data(value, arg);
                    }
                    A & data;
                    value_type const & value;
                };
            public:
                void operator()(value_type const & value, detail::weight_variant_type const & weight) {
                    boost::apply_visitor(call_2_visitor<A>(this->m_data, value), weight);
                }

                void save(hdf5::archive & ar) const { 
                    ar[""] = this->m_data; 
                   }
                void load(hdf5::archive & ar) { 
                    ar[""] >> this->m_data; 
                }

                void print(std::ostream & os) const {
                    this->m_data.print(os);
                }

                void reset() {
                    this->m_data.reset();
                }

#ifdef ALPS_HAVE_MPI
                void collective_merge(
                      boost::mpi::communicator const & comm
                    , int root = 0
                ) {
                    this->m_data.collective_merge(comm, root);
                }

                void collective_merge(
                      boost::mpi::communicator const & comm
                    , int root = 0
                ) const {
                    this->m_data.collective_merge(comm, root);
                }
#endif
        };

        template<typename A> class derived_result_wrapper : public derived_wrapper<A> {
            public:
                derived_result_wrapper(): derived_wrapper<A>() {}

                derived_result_wrapper(A const & arg): derived_wrapper<A>(arg) {}

                base_wrapper<typename value_type<A>::type> * clone() const { 
                    return new derived_result_wrapper<A>(this->m_data); 
                }
                base_wrapper<typename value_type<A>::type> * result() const { 
                    throw std::runtime_error(std::string("A result(") + typeid(A).name() + ") cannot be converted to a result" + ALPS_STACKTRACE);
                    return NULL;
                }

                #define OPERATOR_PROXY(AUGOPNAME, AUGOP)                                            \
                    void AUGOPNAME(base_wrapper<typename value_type<A>::type> const & arg) {        \
                        this->m_data AUGOP arg.template extract<A>();                               \
                    }                                                                               \
                    void AUGOPNAME(double arg) {                                                    \
                        this->m_data AUGOP arg;                                                     \
                    }
                OPERATOR_PROXY(operator+=, +=)
                OPERATOR_PROXY(operator-=, -=)
                OPERATOR_PROXY(operator*=, *=)
                OPERATOR_PROXY(operator/=, /=)
                #undef OPERATOR_PROXY

                void negate() {
                    this->m_data.negate();
                }
                void inverse() {
                    this->m_data.inverse();
                }

                #define FUNCTION_PROXY(FUN)                                                         \
                    void FUN () {                                                                   \
                        this->m_data. FUN ();                                                       \
                    }

                FUNCTION_PROXY(sin)
                FUNCTION_PROXY(cos)
                FUNCTION_PROXY(tan)
                FUNCTION_PROXY(sinh)
                FUNCTION_PROXY(cosh)
                FUNCTION_PROXY(tanh)
                FUNCTION_PROXY(asin)
                FUNCTION_PROXY(acos)
                FUNCTION_PROXY(atan)
                FUNCTION_PROXY(abs)
                FUNCTION_PROXY(sqrt)
                FUNCTION_PROXY(log)
                FUNCTION_PROXY(sq)
                FUNCTION_PROXY(cb)
                FUNCTION_PROXY(cbrt)

                #undef FUNCTION_PROXY
        };
        template<typename T, typename A> derived_result_wrapper<A> operator/(T arg, derived_result_wrapper<A> res) {
            return arg * res.inverse();
        }

        template<typename A> class derived_accumulator_wrapper : public derived_wrapper<A> {
            public:
                derived_accumulator_wrapper(): derived_wrapper<A>() {}

                derived_accumulator_wrapper(A const & arg): derived_wrapper<A>(arg) {}
                
                base_wrapper<typename value_type<A>::type> * clone() const {
                    return new derived_accumulator_wrapper<A>(this->m_data); 
                }
                base_wrapper<typename value_type<A>::type> * result() const { 
                    return result_impl<A>();
                }

                void operator+=(base_wrapper<typename value_type<A>::type> const &) {
                    throw std::runtime_error("The Operator += is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator-=(base_wrapper<typename value_type<A>::type> const &) {
                    throw std::runtime_error("The Operator -= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator*=(base_wrapper<typename value_type<A>::type> const &) {
                    throw std::runtime_error("The Operator *= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator/=(base_wrapper<typename value_type<A>::type> const &) {
                    throw std::runtime_error("The Operator /= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }

                void operator+=(double) {
                    throw std::runtime_error("The operator += is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator-=(double) {
                    throw std::runtime_error("The operator -= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator*=(double) {
                    throw std::runtime_error("The operator *= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void operator/=(double) {
                    throw std::runtime_error("The operator /= is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }

                void negate() {
                    throw std::runtime_error("The function negate is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }
                void inverse() {
                    throw std::runtime_error("The function inverse is not implemented for accumulators, only for results" + ALPS_STACKTRACE);
                }

                #define FUNCTION_PROXY(FUN)                                                                                                           \
                    void FUN () {                                                                                                                     \
                        throw std::runtime_error("The Function " #FUN " is not implemented for accumulators, only for results" + ALPS_STACKTRACE);    \
                    }

                FUNCTION_PROXY(sin)
                FUNCTION_PROXY(cos)
                FUNCTION_PROXY(tan)
                FUNCTION_PROXY(sinh)
                FUNCTION_PROXY(cosh)
                FUNCTION_PROXY(tanh)
                FUNCTION_PROXY(asin)
                FUNCTION_PROXY(acos)
                FUNCTION_PROXY(atan)
                FUNCTION_PROXY(abs)
                FUNCTION_PROXY(sqrt)
                FUNCTION_PROXY(log)
                FUNCTION_PROXY(sq)
                FUNCTION_PROXY(cb)
                FUNCTION_PROXY(cbrt)

                #undef FUNCTION_PROXY

            private:

                template<typename T> typename boost::enable_if<typename has_result_type<T>::type, base_wrapper<typename value_type<A>::type> *>::type result_impl() const {
                    return new derived_result_wrapper<typename A::result_type>(this->m_data);
                }
                template<typename T> typename boost::disable_if<typename has_result_type<T>::type, base_wrapper<typename value_type<A>::type> *>::type result_impl() const {
                    throw std::runtime_error(std::string("The type ") + typeid(A).name() + " has no result_type" + ALPS_STACKTRACE);
                    return NULL;
                }

        };
    }
}

 #endif