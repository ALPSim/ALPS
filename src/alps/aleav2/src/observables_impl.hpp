/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

/**
 * What is behind the classes of <alps/alea/observables.hpp>. Not installed.
 */
#pragma once

#include <alps/alea/observables.hpp>

#include <alps/alea/autocorr.hpp>
#include <alps/alea/batch.hpp>
#include <alps/alea/propagation.hpp>
#include <alps/alea/transform.hpp>
#include <alps/alea/variance.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace alps { namespace alea {

namespace observables_detail {

typedef std::function<std::vector<double>(const std::vector<double> &)> vector_function;

inline std::vector<double> to_vector(const column<double> &in)
{
    return std::vector<double>(in.data(), in.data() + in.size());
}

/** A function of the means of several batch results, for the jackknife */
struct function_transformer
    : public transformer<double>
{
    function_transformer(const vector_function &f, size_t in_size)
        : f_(f)
        , in_size_(in_size)
        , out_size_(f(std::vector<double>(in_size, 1.)).size())
    { }

    column<double> operator()(const column<double> &in) const
    {
        const std::vector<double> out = f_(to_vector(in));
        if (out.size() != out_size_)
            throw std::logic_error("A function of results changed its number of components");
        column<double> res(out.size());
        std::copy(out.begin(), out.end(), res.data());
        return res;
    }

    size_t in_size() const { return in_size_; }

    size_t out_size() const { return out_size_; }

private:
    vector_function f_;
    size_t in_size_, out_size_;
};

inline std::string quoted(const std::string &name)
{
    return name.empty() ? std::string("a function of results") : "'" + name + "'";
}

}

using namespace observables_detail;

struct observable_result::impl
{
    typedef std::shared_ptr<const batch_result<double> > batch_pointer;

    explicit impl(const observable &declaration_)
        : declaration(declaration_)
        , size(0)
        , count(0)
        , convergence(CONVERGED)
    { }

    observable declaration;
    size_t size;
    uint64_t count;

    // What was measured. A histogram has `counts`, an observable without bins
    // `simple`, any other `autocorr` and `batch`.
    std::vector<uint64_t> counts;
    std::shared_ptr<const var_result<double> > simple;
    std::shared_ptr<const autocorr_result<double> > autocorr;
    batch_pointer batch;

    // How mean and error follow from bins: `function` of the means of `parts`,
    // one after the other. An observable with bins is its only part and has no
    // function. A signed observable is the ratio of two parts. A function of
    // results has the parts of its operands, so that it is always evaluated
    // from measured bins and never from the bins of another function.
    std::vector<batch_pointer> parts;
    vector_function function;
    error_convergence convergence;      // worst of the parts, without the own
    mutable batch_pointer evaluated;

    // The binning analysis of `autocorr`: one row per bin length with at
    // least two bins, one column per component. Evaluated once, on demand.
    mutable std::shared_ptr<const std::vector<std::vector<double> > > binning;

    // A result read from a file has the numbers that were stored, and bins if
    // the file has them. It has none of the accumulated data above.
    struct stored_values
    {
        std::vector<double> mean, error, variance, tau;
        std::vector<error_convergence> convergence;
        std::vector<std::vector<double> > binning_errors;    // per bin length
    };
    std::shared_ptr<const stored_values> stored;

    // The samples of a signed observable as measured, sign times quantity,
    // when read from a file. Otherwise they are in `autocorr` and `batch`.
    std::shared_ptr<const impl> raw;

    bool measured() const { return autocorr || simple; }

    void need_samples() const
    {
        if (count == 0)
            throw std::runtime_error("No measurements for " + quoted(declaration.name()));
    }

    void need_values() const
    {
        need_samples();
        if (declaration.histogram_)
            throw std::logic_error(quoted(declaration.name()) + " is a histogram; use histogram()");
    }

    void need_bins() const
    {
        need_values();
        if (parts.empty()) {
            throw std::invalid_argument(
                quoted(declaration.name()) + " has no bins and cannot enter a function of results");
        }
    }

    size_t parts_size() const
    {
        size_t sum = 0;
        for (const batch_pointer &part : parts)
            sum += part->size();
        return sum;
    }

    /** The jackknife over the parts; throws if their bins do not line up */
    const batch_result<double> &evaluate() const
    {
        if (evaluated)
            return *evaluated;

        const batch_data<double> &first = parts.front()->store();
        batch_data<double> joined(parts_size(), first.num_batches());
        joined.count() = first.count();
        joined.set_batch_size(first.batch_size());
        size_t row = 0;
        for (const batch_pointer &part : parts) {
            if (part->store().num_batches() != first.num_batches() || part->store().count() != first.count()) {
                throw std::invalid_argument(
                    "The bins needed for " + quoted(declaration.name()) + " do not line up: its operands"
                    " were not measured the same number of times, or with different numbers of bins");
            }
            joined.batch().middleRows(row, part->size()) = part->store().batch();
            row += part->size();
        }
        const function_transformer tf(function, joined.size());
        evaluated.reset(new batch_result<double>(transform(jackknife_prop(), tf, batch_result<double>(joined))));
        return *evaluated;
    }
};

struct accumulator::impl
{
    explicit impl(const observable &declaration_) : declaration(declaration_), size(0) { reset(); }

    impl(const impl &other)
        : declaration(other.declaration)
        , size(other.size)
        , counts(other.counts)
        , simple(other.simple ? new var_acc<double>(*other.simple) : nullptr)
        , autocorr(other.autocorr ? new autocorr_acc<double>(*other.autocorr) : nullptr)
        , batch(other.batch ? new batch_acc<double>(*other.batch) : nullptr)
    { }

    observable declaration;
    size_t size;                    // 0 until the first sample
    std::vector<uint64_t> counts;
    std::unique_ptr<var_acc<double> > simple;
    std::unique_ptr<autocorr_acc<double> > autocorr;
    std::unique_ptr<batch_acc<double> > batch;

    void reset()
    {
        size = 0;
        simple.reset();
        autocorr.reset();
        batch.reset();
        counts.clear();
        if (declaration.histogram_)
            counts.assign(declaration.max_ - declaration.min_ + 1, 0);
    }

    /** Fix the number of components, at the first sample or the first merge */
    void set_size(size_t new_size)
    {
        if (size == new_size)
            return;
        if (size != 0 || new_size == 0) {
            throw std::invalid_argument(
                quoted(declaration.name()) + " has " + std::to_string(size) + " components, not "
                + std::to_string(new_size));
        }
        size = new_size;
        if (declaration.with_bins_) {
            autocorr.reset(new autocorr_acc<double>(size));
            batch.reset(new batch_acc<double>(size, declaration.num_bins_));
        } else {
            simple.reset(new var_acc<double>(size));
        }
    }

    uint64_t count() const
    {
        if (declaration.histogram_)
            return std::accumulate(counts.begin(), counts.end(), uint64_t(0));
        return autocorr ? autocorr->count() : simple ? simple->count() : 0;
    }

    std::shared_ptr<observable_result::impl> make_result() const
    {
        std::shared_ptr<observable_result::impl> out(new observable_result::impl(declaration));
        out->size = declaration.histogram_ ? counts.size() : size;
        out->count = count();
        out->counts = counts;
        if (simple)
            out->simple.reset(new var_result<double>(simple->result()));
        if (autocorr) {
            out->autocorr.reset(new autocorr_result<double>(autocorr->result()));
            out->batch.reset(new batch_result<double>(batch->result()));
            out->parts.push_back(out->batch);
        }
        return out;
    }
};

struct observable_set::impl
{
    std::map<std::string, accumulator> observables;
};

struct result_set::impl
{
    std::map<std::string, observable_result> results;

    /** Give signed observables the bins of their sign */
    void link(std::map<std::string, std::shared_ptr<observable_result::impl> > &parts)
    {
        results.clear();
        for (auto &named : parts) {
            observable_result::impl &self = *named.second;
            const std::string &sign = self.declaration.sign_;
            if (!sign.empty() && self.batch) {
                if (!parts.count(sign) || !parts[sign]->batch) {
                    throw std::runtime_error(
                        quoted(named.first) + " is signed by " + quoted(sign) + ", which "
                        + (parts.count(sign) ? "has no bins or no measurements" : "is not an observable of the set"));
                }
                const observable_result::impl &sign_part = *parts[sign];
                if (sign_part.size != 1)
                    throw std::runtime_error("The sign " + quoted(sign) + " must have one component");
                self.parts.assign(1, self.batch);
                self.parts.push_back(sign_part.batch);
                const size_t size = self.size;
                self.function = [size](const std::vector<double> &means) {
                    std::vector<double> ratio(means.begin(), means.begin() + size);
                    for (double &value : ratio)
                        value /= means[size];
                    return ratio;
                };
                self.convergence = sign_part.autocorr->converged()[0];
            }
        }
        for (auto &named : parts)
            results[named.first].impl_ = named.second;
    }
};

/** Lets the functions that read and write sets reach behind the classes */
struct observables_access
{
    static std::map<std::string, accumulator> &observables(observable_set &set) { return set.impl_->observables; }
    static const std::map<std::string, accumulator> &observables(const observable_set &set) { return set.impl_->observables; }
    static accumulator::impl &get(accumulator &acc) { return *acc.impl_; }
    static const accumulator::impl &get(const accumulator &acc) { return *acc.impl_; }
    static accumulator make_accumulator(const observable &declaration) { return accumulator(declaration); }
    static const observable_result::impl &get(const observable_result &result) { return *result.impl_; }
    static const std::map<std::string, observable_result> &results(const result_set &set) { return set.impl_->results; }
    static void link(result_set &set, std::map<std::string, std::shared_ptr<observable_result::impl> > &parts)
    {
        set.impl_->link(parts);
    }
    static std::map<std::string, observable_result> &results(result_set &set) { return set.impl_->results; }
    static observable_result make_result(const std::shared_ptr<const observable_result::impl> &impl)
    {
        observable_result result;
        result.impl_ = impl;
        return result;
    }
    static std::string &sign(observable &declaration) { return declaration.sign_; }
    static const std::string &sign(const observable &declaration) { return declaration.sign_; }
    static std::vector<std::string> &labels(observable &declaration) { return declaration.labels_; }
    static const std::vector<std::string> &labels(const observable &declaration) { return declaration.labels_; }
    static std::size_t &num_bins(observable &declaration) { return declaration.num_bins_; }
    static const std::size_t &num_bins(const observable &declaration) { return declaration.num_bins_; }
    static bool &with_bins(observable &declaration) { return declaration.with_bins_; }
    static const bool &with_bins(const observable &declaration) { return declaration.with_bins_; }
    static bool &is_histogram(observable &declaration) { return declaration.histogram_; }
    static const bool &is_histogram(const observable &declaration) { return declaration.histogram_; }
    static std::int64_t &min(observable &declaration) { return declaration.min_; }
    static const std::int64_t &min(const observable &declaration) { return declaration.min_; }
    static std::int64_t &max(observable &declaration) { return declaration.max_; }
    static const std::int64_t &max(const observable &declaration) { return declaration.max_; }
};

}} /* namespace alps::alea */
