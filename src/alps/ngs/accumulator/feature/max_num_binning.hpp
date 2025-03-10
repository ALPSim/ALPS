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

#ifndef ALPS_NGS_ACCUMULATOR_MAX_NUM_BINNING_HPP
#define ALPS_NGS_ACCUMULATOR_MAX_NUM_BINNING_HPP

#include <alps/ngs/accumulator/feature.hpp>
#include <alps/ngs/accumulator/parameter.hpp>
#include <alps/ngs/accumulator/feature/mean.hpp>
#include <alps/ngs/accumulator/feature/count.hpp>
#include <alps/ngs/accumulator/feature/error.hpp>

#include <alps/ngs/numeric.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/ngs/short_print.hpp>

#include <alps/numeric/functional.hpp>
#include <alps/numeric/outer_product.hpp>
#include <alps/type_traits/covariance_type.hpp>

#include <boost/mpl/if.hpp>
#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/type_traits/is_scalar.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <stdexcept>

namespace alps {
    namespace accumulator {
        // this should be called namespace tag { struct max_num_binning; }
        // but gcc <= 4.4 has lookup error, so name it different
        struct max_num_binning_tag;

        namespace detail {
            template<typename C, typename M> class max_num_binning_proxy {
                typedef typename std::size_t size_type;

                public:
                    max_num_binning_proxy(std::vector<M> const & bins, C const & num_elements, size_type const & max_number)
                        : m_max_number(max_number)
                        , m_num_elements(num_elements)
                        , m_bins(bins)
                    {}

                    std::vector<M> const & bins() const {
                        return m_bins;
                    }

                    C num_elements() const {
                        return m_num_elements;
                    }

                    size_type max_number() const {
                        return m_max_number;
                    }

                private:

                    size_type m_max_number;
                    C m_num_elements;
                    std::vector<M> const & m_bins;
            };

            template<typename C, typename M> inline std::ostream & operator<<(std::ostream & os, max_num_binning_proxy<C, M> const & arg) {
                if (arg.bins().empty())
                    os << "No Bins";
                else
                    os << short_print(arg.bins(), 4) << "#" << arg.num_elements();
                return os;
            };
        }

        template<typename T> struct max_num_binning_type {
            typedef detail::max_num_binning_proxy<typename count_type<T>::type, typename mean_type<T>::type> type;
        };

        template<typename T> struct has_feature<T, max_num_binning_tag> {
            template<typename R, typename C> static char helper(R(C::*)() const);
            template<typename C> static char check(boost::integral_constant<std::size_t, sizeof(helper(&C::max_num_binning))>*);
            template<typename C> static double check(...);
            typedef boost::integral_constant<bool, sizeof(char) == sizeof(check<T>(0))> type;
        };

        template<typename T> typename max_num_binning_type<T>::type max_num_binning(T const & arg) {
            return arg.max_num_binning();
        }

        template<typename T> struct covariance_type
            : public alps::covariance_type<typename value_type<T>::type>
        {};

        namespace detail {

            template<typename A> typename boost::enable_if<
                  typename has_feature<A, max_num_binning_tag>::type
                , typename max_num_binning_type<A>::type
            >::type max_num_binning_impl(A const & acc) {
                return max_num_binning(acc);
            }
            template<typename A> typename boost::disable_if<
                  typename has_feature<A, max_num_binning_tag>::type
                , typename max_num_binning_type<A>::type
            >::type max_num_binning_impl(A const & acc) {
                throw std::runtime_error(std::string(typeid(A).name()) + " has no max_num_binning-method" + ALPS_STACKTRACE);
                return *static_cast<typename max_num_binning_type<A>::type *>(NULL);
            }

            template<typename A, typename OP> void transform_impl(
                  A & acc
                , OP op
                , typename boost::enable_if<typename has_feature<A, max_num_binning_tag>::type, int>::type = 0
            ) {
                acc.transform(op);
            }
            template<typename A, typename OP> void transform_impl(
                  A & acc
                , OP op
                , typename boost::disable_if<typename has_feature<A, max_num_binning_tag>::type, int>::type = 0
            ) {
                throw std::runtime_error(std::string(typeid(A).name()) + " has no transform-method" + ALPS_STACKTRACE);
            }
        }

        namespace impl {

            template<typename T, typename B> struct Accumulator<T, max_num_binning_tag, B> : public B {

                public:
                    typedef typename alps::accumulator::max_num_binning_type<B>::type max_num_binning_type;
                    typedef Result<T, max_num_binning_tag, typename B::result_type> result_type;

                    Accumulator()
                        : B()
                        , m_mn_max_number(128)
                        , m_mn_elements_in_bin(0)
                        , m_mn_elements_in_partial(0)
                        , m_mn_partial(T())
                    {}

                    Accumulator(Accumulator const & arg)
                        : B(arg)
                        , m_mn_max_number(arg.m_mn_max_number)
                        , m_mn_elements_in_bin(arg.m_mn_elements_in_bin)
                        , m_mn_elements_in_partial(arg.m_mn_elements_in_partial)
                        , m_mn_partial(arg.m_mn_partial)
                        , m_mn_bins(arg.m_mn_bins)
                    {}

                    template<typename ArgumentPack> Accumulator(ArgumentPack const & args, typename boost::disable_if<is_accumulator<ArgumentPack>, int>::type = 0)
                        : B(args)
                        , m_mn_max_number(args[max_bin_number | 128])
                        , m_mn_elements_in_bin(0)
                        , m_mn_elements_in_partial(0)
                        , m_mn_partial(T())
                    {}

                    max_num_binning_type const max_num_binning() const {
                        return max_num_binning_type(m_mn_bins, m_mn_elements_in_bin, m_mn_max_number);
                    }

                    template <typename OP> void transform(OP) {
                        throw std::runtime_error("Transform can only be applied to a result" + ALPS_STACKTRACE);
                    }

                    template <typename U, typename OP> void transform(U const &, OP) {
                        throw std::runtime_error("Transform can only be applied to a result" + ALPS_STACKTRACE);
                    }

                    using B::operator();
                    void operator()(T const & val) {
                        using alps::ngs::numeric::operator+=;
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator/;
                        using alps::ngs::numeric::detail::check_size;

                        B::operator()(val);

                        if (!m_mn_elements_in_bin) {
                            m_mn_bins.push_back(val);
                            m_mn_elements_in_bin = 1;
                        } else {
                            check_size(m_mn_bins[0], val);
                            check_size(m_mn_partial, val);
                            m_mn_partial += val;
                            ++m_mn_elements_in_partial;
                        }

                        // TODO: make library for scalar type
                        typename alps::hdf5::scalar_type<T>::type elements_in_bin = m_mn_elements_in_bin;
                        typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type two = 2;

                        if (m_mn_elements_in_partial == m_mn_elements_in_bin && m_mn_bins.size() >= m_mn_max_number) {
                            if (m_mn_max_number % 2 == 1) {
                                m_mn_partial += m_mn_bins[m_mn_max_number - 1];
                                m_mn_elements_in_partial += m_mn_elements_in_bin;
                            }
                            for (typename count_type<T>::type i = 0; i < m_mn_max_number / 2; ++i)
                                m_mn_bins[i] = (m_mn_bins[2 * i] + m_mn_bins[2 * i + 1]) / two;
                            m_mn_bins.erase(m_mn_bins.begin() + m_mn_max_number / 2, m_mn_bins.end());
                            m_mn_elements_in_bin *= (typename count_type<T>::type)2;
                        }
                        if (m_mn_elements_in_partial == m_mn_elements_in_bin) {
                            m_mn_bins.push_back(m_mn_partial / elements_in_bin);
                            m_mn_partial = T();
                            m_mn_elements_in_partial = 0;
                        }
                    }

                    template<typename S> void print(S & os) const {
                        B::print(os);
                        os << " Bins: " << max_num_binning();
                    }

                    void save(hdf5::archive & ar) const {
                        B::save(ar);
                        if (B::count()) {
                            ar["timeseries/partialbin"] = m_mn_partial;
                            ar["timeseries/partialbin/@count"] = m_mn_elements_in_partial;
                        }
                        ar["timeseries/data"] = m_mn_bins;
                        ar["timeseries/data/@binningtype"] = "linear";
                        ar["timeseries/data/@minbinsize"] = 0; // TODO: what should we put here?
                        ar["timeseries/data/@binsize"] = m_mn_elements_in_bin;
                        ar["timeseries/data/@maxbinnum"] = m_mn_max_number;
                    }

                    void load(hdf5::archive & ar) { // TODO: make archive const
                        B::load(ar);
                        ar["timeseries/data"] >> m_mn_bins;
                        ar["timeseries/data/@binsize"] >> m_mn_elements_in_bin;
                        ar["timeseries/data/@maxbinnum"] >> m_mn_max_number;
                        if (ar.is_data("timeseries/partialbin")) {
                            ar["timeseries/partialbin"] >> m_mn_partial;
                            ar["timeseries/partialbin/@count"] >> m_mn_elements_in_partial;
                        }
                    }

                    static std::size_t rank() { return B::rank() + 1; }
                    static bool can_load(hdf5::archive & ar) { // TODO: make archive const
                        using alps::hdf5::get_extent;

                        return B::can_load(ar)
                            && ar.is_data("timeseries/data")
                            && ar.is_attribute("timeseries/data/@binsize")
                            && ar.is_attribute("timeseries/data/@maxbinnum")
                            && get_extent(T()).size() + 1 == ar.dimensions("timeseries/data")
                        ;
                    }

                    void reset() {
                        B::reset();
                        m_mn_elements_in_bin = typename B::count_type();
                        m_mn_elements_in_partial = typename B::count_type();
                        m_mn_partial = T();
                        m_mn_bins = std::vector<typename mean_type<B>::type>();
                    }

#ifdef ALPS_HAVE_MPI
                    void collective_merge(
                          boost::mpi::communicator const & comm
                        , int root
                    ) {
                        if (comm.rank() == root) {
                            B::collective_merge(comm, root);
                            if (!m_mn_bins.empty()) {
                                std::vector<typename mean_type<B>::type> local_bins(m_mn_bins), merged_bins;
                                partition_bins(comm, local_bins, merged_bins, root);
                                B::reduce_if(comm, merged_bins, m_mn_bins, std::plus<typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type>(), root);
                            }
                        } else
                            const_cast<Accumulator<T, max_num_binning_tag, B> const *>(this)->collective_merge(comm, root);
                    }

                    void collective_merge(
                          boost::mpi::communicator const & comm
                        , int root
                    ) const {
                        B::collective_merge(comm, root);
                        if (comm.rank() == root)
                            throw std::runtime_error("A const object cannot be root" + ALPS_STACKTRACE);
                        else if (!m_mn_bins.empty()) {
                            std::vector<typename mean_type<B>::type> local_bins(m_mn_bins), merged_bins;
                            partition_bins(comm, local_bins, merged_bins, root);
                            B::reduce_if(comm, merged_bins, std::plus<typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type>(), root);
                        }
                    }

                private:
                    void partition_bins (
                          boost::mpi::communicator const & comm
                        , std::vector<typename mean_type<B>::type> & local_bins
                        , std::vector<typename mean_type<B>::type> & merged_bins
                        , int root
                    ) const {
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator/;
                        using alps::ngs::numeric::detail::check_size;

                        typename B::count_type elements_in_local_bins = boost::mpi::all_reduce(comm, m_mn_elements_in_bin, boost::mpi::maximum<typename B::count_type>());
                        typename B::count_type howmany = (elements_in_local_bins - 1) / m_mn_elements_in_bin + 1;
                        if (howmany > 1) {
                            typename B::count_type newbins = local_bins.size() / howmany;
                            typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type howmany_vt = howmany;
                            for (typename B::count_type i = 0; i < newbins; ++i) {
                                local_bins[i] = local_bins[howmany * i];
                                for (typename B::count_type j = 1; j < howmany; ++j)
                                    local_bins[i] = local_bins[i] + local_bins[howmany * i + j];
                                local_bins[i] = local_bins[i] / howmany_vt;
                            }
                            local_bins.resize(newbins);
                        }

                        std::vector<std::size_t> index(comm.size());
                        boost::mpi::all_gather(comm, local_bins.size(), index);
                        std::size_t total_bins = std::accumulate(index.begin(), index.end(), 0);
                        std::size_t perbin = total_bins < m_mn_max_number ? 1 : total_bins / m_mn_max_number;
                        typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type perbin_vt = perbin;

                        merged_bins.resize(perbin == 1 ? total_bins : m_mn_max_number);
                        for (typename std::vector<typename mean_type<B>::type>::iterator it = merged_bins.begin(); it != merged_bins.end(); ++it)
                            check_size(*it, local_bins[0]);

                        std::size_t start = std::accumulate(index.begin(), index.begin() + comm.rank(), 0);
                        for (std::size_t i = start / perbin, j = start % perbin, k = 0; i < merged_bins.size() && k < local_bins.size(); ++k) {
                            merged_bins[i] = merged_bins[i] + local_bins[k] / perbin_vt;
                            if (++j == perbin)
                                ++i, j = 0;
                        }
                    }
#endif

                private:

                    std::size_t m_mn_max_number;
                    typename B::count_type m_mn_elements_in_bin, m_mn_elements_in_partial;
                    T m_mn_partial;
                    std::vector<typename mean_type<B>::type> m_mn_bins;
            };

            template<typename T, typename B> class Result<T, max_num_binning_tag, B> : public B {

                public:
                    typedef typename alps::accumulator::max_num_binning_type<B>::type max_num_binning_type;

                    Result()
                        : B()
                        , m_mn_max_number(0) 
                        , m_mn_elements_in_bin(0)
                        , m_mn_count(typename B::count_type())
                        , m_mn_mean(typename mean_type<B>::type())
                        , m_mn_error(typename error_type<B>::type())
                        , m_mn_cannot_rebin(false)
                        , m_mn_jackknife_valid(false)
                        , m_mn_data_is_analyzed(true)
                        , m_mn_jackknife_bins(0)
                    {}

                    // copy constructor
                    template<typename A> Result(A const & acc, typename boost::enable_if<boost::is_base_of<ResultBase<T>, A>, int>::type = 0)
                        : B(acc)
                        , m_mn_max_number(acc.m_mn_max_number)
                        , m_mn_elements_in_bin(acc.m_mn_elements_in_bin)
                        , m_mn_bins(acc.m_mn_bins)
                        , m_mn_count(acc.count())
                        , m_mn_mean(acc.mean())
                        , m_mn_error(acc.error())
                        , m_mn_cannot_rebin(acc.m_mn_cannot_rebin)
                        , m_mn_jackknife_valid(acc.m_mn_jackknife_valid)
                        , m_mn_data_is_analyzed(acc.m_mn_data_is_analyzed)
                        , m_mn_jackknife_bins(acc.m_mn_jackknife_bins)
                    {}

                    template<typename A> Result(A const & acc, typename boost::disable_if<boost::is_base_of<ResultBase<T>, A>, int>::type = 0)
                        : B(acc)
                        , m_mn_max_number(detail::max_num_binning_impl(acc).max_number())
                        , m_mn_elements_in_bin(detail::max_num_binning_impl(acc).num_elements())
                        , m_mn_bins(detail::max_num_binning_impl(acc).bins())
                        , m_mn_count(acc.count())
                        , m_mn_mean(acc.mean())
                        , m_mn_error(acc.error())
                        , m_mn_cannot_rebin(false)
                        , m_mn_jackknife_valid(false)
                        , m_mn_data_is_analyzed(true)
                        , m_mn_jackknife_bins(0)
                    {}

                    typename B::count_type count() const { 
                        analyze();
                        return m_mn_count;
                    }

                    typename mean_type<B>::type const & mean() const {
                        analyze();
                        return m_mn_mean;
                    }

                    typename error_type<B>::type const & error() const {
                        analyze();
                        return m_mn_error;
                    }

                    max_num_binning_type const max_num_binning() const {
                        return max_num_binning_type(m_mn_bins, m_mn_elements_in_bin, m_mn_max_number);
                    }

                    template <typename A> typename boost::enable_if<
                        typename has_feature<A, max_num_binning_tag>::type, typename covariance_type<B>::type
                    >::type covariance(A const & obs) const {
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator/;
                        using alps::numeric::outer_product;

                        generate_jackknife();
                        obs.generate_jackknife();
                        if (m_mn_jackknife_bins.size() != obs.m_mn_jackknife_bins.size())
                            throw std::runtime_error("Unequal number of bins in calculation of covariance matrix" + ALPS_STACKTRACE);
                        if (!m_mn_jackknife_bins.size() || !obs.m_mn_jackknife_bins.size())
                            throw std::runtime_error("No binning information available for calculation of covariances" + ALPS_STACKTRACE);

                        typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type bin_number = m_mn_bins.size();

                        typename mean_type<B>::type unbiased_mean_1;
                        for (typename std::vector<typename mean_type<B>::type>::const_iterator it = m_mn_jackknife_bins.begin() + 1; it != m_mn_jackknife_bins.end(); ++it)
                            unbiased_mean_1 = unbiased_mean_1 + *it / bin_number;

                        typename mean_type<B>::type unbiased_mean_2;
                        for (typename std::vector<typename mean_type<B>::type>::const_iterator it = obs.m_mn_jackknife_bins.begin() + 1; it != obs.m_mn_jackknife_bins.end(); ++it)
                            unbiased_mean_2 = unbiased_mean_2 + *it / bin_number;

                        typename covariance_type<B>::type cov = outer_product(m_mn_jackknife_bins[1], obs.m_mn_jackknife_bins[1]);
                        for (typename B::count_type i = 1; i < m_mn_bins.size(); ++i)
                            cov += outer_product(m_mn_jackknife_bins[i + 1], obs.m_mn_jackknife_bins[i + 1]);
                        cov /= bin_number;
                        cov -= outer_product(unbiased_mean_1, unbiased_mean_2);
                        cov *= bin_number - 1;
                        return cov;
                    }

                    // Adapted from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance
                    // It is a two-pass algorith, which first calculates estimates for the mean and then performs
                    // the stable algorithm on the residuals. According to literature and local authorities, this
                    // is the most accurate and stable way to calculate variances.
                    template <typename A> typename boost::enable_if<
                        typename has_feature<A, max_num_binning_tag>::type, typename covariance_type<B>::type
                    >::type accurate_covariance(A const & obs) const {
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator-;
                        using alps::ngs::numeric::operator/;
                        using alps::numeric::outer_product;

                        generate_jackknife();
                        obs.generate_jackknife();
                        if (m_mn_jackknife_bins.size() != obs.m_mn_jackknife_bins.size())
                            throw std::runtime_error("Unequal number of bins in calculation of covariance matrix" + ALPS_STACKTRACE);
                        if (!m_mn_jackknife_bins.size() || !obs.m_mn_jackknife_bins.size())
                            throw std::runtime_error("No binning information available for calculation of covariances" + ALPS_STACKTRACE);

                        typedef typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type scalar_type;
                        scalar_type bin_number = m_mn_bins.size();

                        typename mean_type<B>::type unbiased_mean_1;
                        for (typename std::vector<typename mean_type<B>::type>::const_iterator it = m_mn_jackknife_bins.begin() + 1; it != m_mn_jackknife_bins.end(); ++it)
                            unbiased_mean_1 = unbiased_mean_1 + *it / bin_number;

                        typename mean_type<B>::type unbiased_mean_2;
                        for (typename std::vector<typename mean_type<B>::type>::const_iterator it = obs.m_mn_jackknife_bins.begin() + 1; it != obs.m_mn_jackknife_bins.end(); ++it)
                            unbiased_mean_2 = unbiased_mean_2 + *it / bin_number;

                        std::vector<typename mean_type<B>::type> X(m_mn_bins.size());
                        std::vector<typename mean_type<B>::type> Y(m_mn_bins.size());
                        for (typename B::count_type i = 0; i < m_mn_bins.size(); ++i) {
                            X[i] = m_mn_jackknife_bins[i + 1] - unbiased_mean_1;
                            Y[i] = obs.m_mn_jackknife_bins[i + 1] - unbiased_mean_2;
                        }

                        typename mean_type<B>::type xbar;
                        typename mean_type<B>::type ybar;
                        typename covariance_type<B>::type cov = outer_product(xbar, ybar);
                        for (typename B::count_type i = 0; i < m_mn_bins.size(); ++i) {
                            typename mean_type<B>::type delta_x = X[i] - xbar;
                            typename mean_type<B>::type delta_y = Y[i] - ybar;
                            xbar = xbar + delta_x / scalar_type(i + 1);
                            cov += outer_product(X[i] - xbar, delta_y);
                            ybar = ybar + delta_y / scalar_type(i + 1);
                        }
                        cov /= bin_number;
                        cov *= bin_number - 1;
                        return cov;
                    }

                    // TODO: use mean error from here ...
                    template<typename S> void print(S & os) const {
                        // TODO: use m_mn_variables!
                        B::print(os);
                        os << " Bins: " << max_num_binning();
                    }

                    void save(hdf5::archive & ar) const {
                        analyze();
                        ar["count"] = m_mn_count;
                        ar["@cannotrebin"] = m_mn_cannot_rebin;
                        ar["mean/value"] = m_mn_mean;
                        ar["mean/error"] = m_mn_error;
                        ar["timeseries/data"] = m_mn_bins;
                        ar["timeseries/data/@binsize"] = m_mn_elements_in_bin;
                        ar["timeseries/data/@maxbinnum"] = m_mn_max_number;
                        ar["timeseries/data/@binningtype"] = "linear";
                        if (m_mn_jackknife_valid) {
                            ar["jacknife/data"] = m_mn_jackknife_bins;
                            ar["jacknife/data/@binningtype"] = "linear";
                        }
                    }

                    void load(hdf5::archive & ar) {
                        // TODO: implement imformations above
                        B::load(ar);
                        ar["timeseries/data"] >> m_mn_bins;
                        ar["timeseries/data/@binsize"] >> m_mn_elements_in_bin;
                        ar["timeseries/data/@maxbinnum"] >> m_mn_max_number;
                    }

                    static std::size_t rank() { return B::rank() + 1; }
                    static bool can_load(hdf5::archive & ar) { // TODO: make archive const
                        using alps::hdf5::get_extent;

                        return B::can_load(ar)
                            && ar.is_data("timeseries/data")
                            && ar.is_attribute("timeseries/data/@binsize")
                            && ar.is_attribute("timeseries/data/@maxbinnum")
                            && get_extent(T()).size() + 1 == ar.dimensions("timeseries/data")
                        ;
                    }

                    template<typename U> void operator+=(U const & arg) { augadd(arg); }
                    template<typename U> void operator-=(U const & arg) { augsub(arg); }
                    template<typename U> void operator*=(U const & arg) { augmul(arg); }
                    template<typename U> void operator/=(U const & arg) { augdiv(arg); }

                    template <typename OP> void transform(OP op) {
                        generate_jackknife();
                        m_mn_data_is_analyzed = false;
                        m_mn_cannot_rebin = true;
                        typename std::vector<typename mean_type<B>::type>::iterator it;
                        for (it = m_mn_bins.begin(); it != m_mn_bins.end(); ++it)
                            *it = op(*it);
                        for (it = m_mn_jackknife_bins.begin(); it != m_mn_jackknife_bins.end(); ++it)
                            *it = op(*it);
                        analyze();
                    }

                    template <typename OP, typename U> void transform(OP op, U const & arg) {
                        generate_jackknife();
                        arg.generate_jackknife(); /* TODO: make this more generic */
                        if (arg.m_mn_jackknife_bins.size() != m_mn_jackknife_bins.size()) /* TODO: make this more generic */
                            throw std::runtime_error("Unable to transform: unequal number of bins" + ALPS_STACKTRACE);
                        m_mn_data_is_analyzed = false;
                        m_mn_cannot_rebin = true;
                        typename std::vector<typename mean_type<B>::type>::iterator it;
                        typename std::vector<typename mean_type<U>::type>::const_iterator jt;
                        for (it = m_mn_bins.begin(), jt = arg.m_mn_bins.begin(); it != m_mn_bins.end(); ++it, ++jt)
                            *it = op(*it, *jt);
                        for (it = m_mn_jackknife_bins.begin(), jt = arg.m_mn_jackknife_bins.begin(); it != m_mn_jackknife_bins.end(); ++it, ++jt)
                            *it = op(*it, *jt);
                    }

                    #define NUMERIC_FUNCTION_IMPLEMENTATION(FUNCTION_NAME)                                                              \
                        void FUNCTION_NAME () {                                                                                         \
                            using alps::numeric::sq;                                                                                    \
                            using alps::ngs::numeric::sq;                                                                               \
                            using alps::numeric::cbrt;                                                                                  \
                            using alps::ngs::numeric::cbrt;                                                                             \
                            using alps::numeric::cb;                                                                                    \
                            using alps::ngs::numeric::cb;                                                                               \
                            using std::sqrt;                                                                                            \
                            using alps::ngs::numeric::sqrt;                                                                             \
                            using std::exp;                                                                                             \
                            using alps::ngs::numeric::exp;                                                                              \
                            using std::log;                                                                                             \
                            using alps::ngs::numeric::log;                                                                              \
                            using std::abs;                                                                                             \
                            using alps::ngs::numeric::abs;                                                                              \
                            using std::pow;                                                                                             \
                            using alps::ngs::numeric::pow;                                                                              \
                            using std::sin;                                                                                             \
                            using alps::ngs::numeric::sin;                                                                              \
                            using std::cos;                                                                                             \
                            using alps::ngs::numeric::cos;                                                                              \
                            using std::tan;                                                                                             \
                            using alps::ngs::numeric::tan;                                                                              \
                            using std::sinh;                                                                                            \
                            using alps::ngs::numeric::sinh;                                                                             \
                            using std::cosh;                                                                                            \
                            using alps::ngs::numeric::cosh;                                                                             \
                            using std::tanh;                                                                                            \
                            using alps::ngs::numeric::tanh;                                                                             \
                            using std::asin;                                                                                            \
                            using alps::ngs::numeric::asin;                                                                             \
                            using std::acos;                                                                                            \
                            using alps::ngs::numeric::acos;                                                                             \
                            using std::atan;                                                                                            \
                            using alps::ngs::numeric::atan;                                                                             \
                            transform((typename value_type<B>::type(*)(typename value_type<B>::type))& FUNCTION_NAME );                 \
                            B:: FUNCTION_NAME ();                                                                                       \
                        }

                    NUMERIC_FUNCTION_IMPLEMENTATION(sin)
                    NUMERIC_FUNCTION_IMPLEMENTATION(cos)
                    NUMERIC_FUNCTION_IMPLEMENTATION(tan)
                    NUMERIC_FUNCTION_IMPLEMENTATION(sinh)
                    NUMERIC_FUNCTION_IMPLEMENTATION(cosh)
                    NUMERIC_FUNCTION_IMPLEMENTATION(tanh)
                    NUMERIC_FUNCTION_IMPLEMENTATION(asin)
                    NUMERIC_FUNCTION_IMPLEMENTATION(acos)
                    NUMERIC_FUNCTION_IMPLEMENTATION(atan)
                    NUMERIC_FUNCTION_IMPLEMENTATION(abs)
                    NUMERIC_FUNCTION_IMPLEMENTATION(sq)
                    NUMERIC_FUNCTION_IMPLEMENTATION(sqrt)
                    NUMERIC_FUNCTION_IMPLEMENTATION(cb)
                    NUMERIC_FUNCTION_IMPLEMENTATION(cbrt)
                    NUMERIC_FUNCTION_IMPLEMENTATION(exp)
                    NUMERIC_FUNCTION_IMPLEMENTATION(log)

                    #undef NUMERIC_FUNCTION_IMPLEMENTATION

                private:
                    std::size_t m_mn_max_number;
                    typename B::count_type m_mn_elements_in_bin;
                    std::vector<typename mean_type<B>::type> m_mn_bins;
                    mutable typename B::count_type m_mn_count;
                    mutable typename mean_type<B>::type m_mn_mean;
                    mutable typename error_type<B>::type m_mn_error;
                    mutable bool m_mn_cannot_rebin;
                    mutable bool m_mn_jackknife_valid;
                    mutable bool m_mn_data_is_analyzed;
                    mutable std::vector<typename mean_type<B>::type> m_mn_jackknife_bins;

                    void generate_jackknife() const {
                        using alps::ngs::numeric::operator-;
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator*;
                        using alps::ngs::numeric::operator/;
                        // build jackknife data structure
                        if (!m_mn_bins.empty() && !m_mn_jackknife_valid) {
                            if (m_mn_cannot_rebin)
                                throw std::runtime_error("Cannot build jackknife data structure after nonlinear operations" + ALPS_STACKTRACE);
                            m_mn_jackknife_bins.clear();
                            m_mn_jackknife_bins.resize(m_mn_bins.size() + 1);
                            // Order-N initialization of jackknife data structure
                            //    m_mn_jackknife_bins[0]   =  <x>
                            //    m_mn_jackknife_bins[i+1] =  <x_i>_{jacknife}
                            typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type bin_number = m_mn_bins.size();
                            for(std::size_t j = 0; j < m_mn_bins.size(); ++j) // m_mn_jackknife_bins[0] = \sum_{j} m_mn_bins[j]
                                m_mn_jackknife_bins[0] = m_mn_jackknife_bins[0] + m_mn_bins[j];
                            for(std::size_t i = 0; i < m_mn_bins.size(); ++i) // m_mn_jackknife_bins[i+1] = \sum_{j != i} m_mn_bins[j] / #m_mn_bins
                                m_mn_jackknife_bins[i + 1] = (m_mn_jackknife_bins[0] - m_mn_bins[i]) / (bin_number - 1);
                            m_mn_jackknife_bins[0] = m_mn_jackknife_bins[0] / bin_number; // m_mn_jackknife_bins[0] is the jacknife mean...
                        }
                        m_mn_jackknife_valid = true;
                    }

                    void analyze() const {
                        using alps::numeric::sq;
                        using std::sqrt;
                        using alps::numeric::sqrt;
                        using alps::ngs::numeric::operator-;
                        using alps::ngs::numeric::operator+;
                        using alps::ngs::numeric::operator*;
                        using alps::ngs::numeric::operator/;
                        if (m_mn_bins.empty())
                            throw std::runtime_error("No Measurement" + ALPS_STACKTRACE);
                        if (!m_mn_data_is_analyzed) {
                            m_mn_count = m_mn_elements_in_bin * m_mn_bins.size();
                            generate_jackknife();
                            if (m_mn_jackknife_bins.size()) {
                                typename mean_type<B>::type unbiased_mean = typename mean_type<B>::type();
                                typename alps::hdf5::scalar_type<typename mean_type<B>::type>::type bin_number = m_mn_bins.size();
                                for (typename std::vector<typename mean_type<B>::type>::const_iterator it = m_mn_jackknife_bins.begin() + 1; it != m_mn_jackknife_bins.end(); ++it)
                                    unbiased_mean = unbiased_mean + *it / bin_number;
                                m_mn_mean = m_mn_jackknife_bins[0] - (unbiased_mean - m_mn_jackknife_bins[0]) * (bin_number - 1);
                                m_mn_error = typename error_type<B>::type();
                                for (std::size_t i = 0; i < m_mn_bins.size(); ++i)
                                    m_mn_error = m_mn_error + sq(m_mn_jackknife_bins[i + 1] - unbiased_mean);
                                m_mn_error = sqrt(m_mn_error / bin_number * (bin_number - 1));
                            }
                        }
                        m_mn_data_is_analyzed = true;
                    }

                    #define NUMERIC_FUNCTION_OPERATOR(OP_NAME, OPEQ_NAME, OP, OP_TOKEN, OP_STD)                                                                                 \
                        template<typename U> void aug ## OP_TOKEN (U const & arg, typename boost::disable_if<boost::is_scalar<U>, int>::type = 0) {                             \
                            typedef typename value_type<B>::type self_value_type;                                                                                               \
                            typedef typename value_type<U>::type arg_value_type;                                                                                                \
                            transform(boost::function<self_value_type(self_value_type, arg_value_type)>( OP_STD <self_value_type, self_value_type, arg_value_type>()), arg);    \
                            B:: OPEQ_NAME (arg);                                                                                                                                \
                        }                                                                                                                                                       \
                        template<typename U> void aug ## OP_TOKEN (U const & arg, typename boost::enable_if<boost::mpl::and_<                                                   \
                              boost::is_scalar<U>                                                                                                                               \
                            , typename has_operator_ ## OP_TOKEN <typename mean_type<B>::type, U>::type                                                                         \
                        >, int>::type = 0) {                                                                                                                                    \
                            using alps::ngs::numeric:: OP_NAME ;                                                                                                                \
                            generate_jackknife();                                                                                                                               \
                            m_mn_data_is_analyzed = false;                                                                                                                      \
                            m_mn_cannot_rebin = true;                                                                                                                           \
                            typename std::vector<typename mean_type<B>::type>::iterator it;                                                                                     \
                            for (it = m_mn_bins.begin(); it != m_mn_bins.end(); ++it)                                                                                           \
                                *it = *it OP arg;                                                                                                                               \
                            for (it = m_mn_jackknife_bins.begin(); it != m_mn_jackknife_bins.end(); ++it)                                                                       \
                                *it = *it OP arg;                                                                                                                               \
                            analyze();                                                                                                                                          \
                            B:: OPEQ_NAME (arg);                                                                                                                                \
                        }                                                                                                                                                       \
                        template<typename U> void aug ## OP_TOKEN (U const & arg, typename boost::enable_if<boost::mpl::and_<                                                   \
                              boost::is_scalar<U>                                                                                                                               \
                            , boost::mpl::not_<typename has_operator_ ## OP_TOKEN <typename mean_type<B>::type, U>::type>                                                       \
                        >, int>::type = 0) {                                                                                                                                    \
                            throw std::runtime_error(std::string(typeid(typename mean_type<B>::type).name())                                                                    \
                                + " has no operator " #OP " with " + typeid(U).name() + ALPS_STACKTRACE);                                                                       \
                        }

                    NUMERIC_FUNCTION_OPERATOR(operator+, operator+=, +, add, alps::numeric::plus)
                    NUMERIC_FUNCTION_OPERATOR(operator-, operator-=, -, sub, alps::numeric::minus)
                    NUMERIC_FUNCTION_OPERATOR(operator*, operator*=, *, mul, alps::numeric::multiplies)
                    NUMERIC_FUNCTION_OPERATOR(operator/, operator/=, /, div, alps::numeric::divides)

                    #undef NUMERIC_FUNCTION_OPERATOR
            };

            template<typename T, typename B> class BaseWrapper<T, max_num_binning_tag, B> : public B {
                public:
                    virtual bool has_max_num_binning() const = 0;
                    virtual typename max_num_binning_type<B>::type max_num_binning() const = 0;
                    virtual void transform(boost::function<typename value_type<B>::type(typename value_type<B>::type)>) = 0;
            };

            template<typename T, typename B> class DerivedWrapper<T, max_num_binning_tag, B> : public B {
                public:
                    DerivedWrapper(): B() {}
                    DerivedWrapper(T const & arg): B(arg) {}

                    bool has_max_num_binning() const { return has_feature<T, max_num_binning_tag>::type::value; }
                    bool has_transform() const { return has_feature<T, max_num_binning_tag>::type::value; }

                    typename max_num_binning_type<B>::type max_num_binning() const { return detail::max_num_binning_impl(this->m_data); }
                    void transform(boost::function<typename value_type<B>::type(typename value_type<B>::type)> op) { return detail::transform_impl(this->m_data, op); }
            };

        }
    }
}

 #endif
