/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include "observables_impl.hpp"

#include <memory>
#include <ostream>
#include <vector>

namespace alps { namespace alea {


// DECLARATION

observable::observable(const std::string &name)
    : name_(name)
    , num_bins_(128)
    , components_(0)
    , with_bins_(true)
    , histogram_(false)
    , vector_(false)
    , min_(0)
    , max_(0)
{ }

observable &observable::components(std::size_t size)
{
    components_ = size;
    return *this;
}

observable &observable::without_bins()
{
    with_bins_ = false;
    return *this;
}

observable &observable::bins(std::size_t num_bins)
{
    if (num_bins == 0 || num_bins % 2 != 0)
        throw std::invalid_argument("The number of bins of " + quoted(name_) + " must be even");
    num_bins_ = num_bins;
    return *this;
}

observable &observable::signed_by(const std::string &sign)
{
    sign_ = sign;
    return *this;
}

observable &observable::labels(const std::vector<std::string> &labels)
{
    labels_ = labels;
    return *this;
}

observable histogram(const std::string &name, std::int64_t min, std::int64_t max)
{
    if (max < min)
        throw std::invalid_argument("The histogram " + quoted(name) + " has no bins");
    observable declaration(name);
    declaration.with_bins_ = false;
    declaration.histogram_ = true;
    declaration.min_ = min;
    declaration.max_ = max;
    return declaration;
}

// RESULT

observable_result::observable_result() : impl_(new impl(observable(""))) { }
observable_result::observable_result(const observable_result &other) : impl_(other.impl_) { }
observable_result &observable_result::operator=(const observable_result &other) { impl_ = other.impl_; return *this; }
observable_result::~observable_result() { }

std::size_t observable_result::size() const { return impl_->size; }

std::uint64_t observable_result::count() const { return impl_->count; }

std::vector<double> observable_result::mean_vector() const
{
    impl_->need_values();
    if (impl_->stored)
        return impl_->stored->mean;
    if (impl_->function)
        return to_vector(impl_->evaluate().mean());
    return to_vector(impl_->autocorr ? impl_->autocorr->mean() : impl_->simple->mean());
}

std::vector<double> observable_result::error_vector() const
{
    impl_->need_values();
    if (impl_->stored)
        return impl_->stored->error;
    if (impl_->function)
        return to_vector(impl_->evaluate().stderror());
    return to_vector(impl_->autocorr ? impl_->autocorr->stderror() : impl_->simple->stderror());
}

std::vector<double> observable_result::variance_vector() const
{
    impl_->need_values();
    if (impl_->stored)
        return impl_->stored->variance;
    if (impl_->function)
        return std::vector<double>(impl_->size, std::numeric_limits<double>::quiet_NaN());
    return to_vector(impl_->autocorr ? impl_->autocorr->level(0).var() : impl_->simple->var());
}

std::vector<double> observable_result::tau_vector() const
{
    impl_->need_values();
    if (impl_->stored)
        return impl_->stored->tau;
    if (impl_->function || !impl_->autocorr)
        return std::vector<double>(impl_->size, std::numeric_limits<double>::quiet_NaN());
    return to_vector(impl_->autocorr->tau());
}

std::vector<error_convergence> observable_result::converged_vector() const
{
    impl_->need_values();
    if (impl_->stored)
        return impl_->stored->convergence;
    std::vector<error_convergence> conv(impl_->size, impl_->convergence);
    if (impl_->autocorr) {
        const std::vector<error_convergence> own = impl_->autocorr->converged();
        for (size_t i = 0; i != conv.size(); ++i)
            conv[i] = std::max(conv[i], own[i]);
    }
    return conv;
}

namespace {

template <typename T>
T only(const std::vector<T> &values, const char *what)
{
    if (values.size() != 1) {
        throw std::logic_error(std::string(what) + "() is for results with one component, this one has "
                               + std::to_string(values.size()));
    }
    return values[0];
}

}

double observable_result::mean() const { return only(mean_vector(), "mean"); }
double observable_result::error() const { return only(error_vector(), "error"); }
double observable_result::tau() const { return only(tau_vector(), "tau"); }
error_convergence observable_result::converged() const { return only(converged_vector(), "converged"); }

const std::vector<std::string> &observable_result::labels() const { return impl_->declaration.labels_; }

bool observable_result::has_bins() const { return !impl_->parts.empty(); }

std::uint64_t observable_result::bin_size() const
{
    impl_->need_bins();
    return impl_->parts.front()->store().batch_size();
}

std::vector<std::vector<double> > observable_result::bin_means() const
{
    impl_->need_bins();
    if (!impl_->batch)
        throw std::logic_error("A function of results has no bins of its own");
    const batch_data<double> &store = impl_->batch->store();
    std::vector<std::vector<double> > bins;
    for (size_t i = 0; i != store.num_batches(); ++i) {
        if (store.count()(i) >= store.batch_size())
            bins.push_back(to_vector(store.batch().col(i) / double(store.count()(i))));
    }
    return bins;
}

std::vector<std::vector<double> > observable_result::covariance() const
{
    const std::vector<std::vector<double> > bins = bin_means();
    const size_t m = bins.size(), n = size();
    if (m < 2)
        throw std::runtime_error(quoted(impl_->declaration.name()) + " has fewer than two full bins");
    std::vector<double> mean(n, 0.);
    for (const std::vector<double> &bin : bins)
        for (size_t k = 0; k != n; ++k)
            mean[k] += bin[k] / m;
    std::vector<std::vector<double> > cov(n, std::vector<double>(n, 0.));
    for (const std::vector<double> &bin : bins)
        for (size_t k = 0; k != n; ++k)
            for (size_t l = 0; l != n; ++l)
                cov[k][l] += (bin[k] - mean[k]) * (bin[l] - mean[l]) / (double(m) * (m - 1));
    return cov;
}

std::vector<double> observable_result::binning_errors(std::size_t i) const
{
    std::vector<double> errors;
    if (impl_->stored) {
        for (const std::vector<double> &level : impl_->stored->binning_errors)
            errors.push_back(level.at(i));
        return errors;
    }
    // A ratio or another function of results has no binning analysis.
    if (!impl_->autocorr || impl_->function)
        return errors;
    // Levels with one batch only say nothing. Each level is evaluated for all
    // components at once, so asking for every component stays linear.
    if (!impl_->binning) {
        std::shared_ptr<std::vector<std::vector<double> > > table(new std::vector<std::vector<double> >);
        for (size_t level = 0; level != impl_->autocorr->nlevel(); ++level) {
            if (impl_->autocorr->level(level).observations() >= 2)
                table->push_back(to_vector(impl_->autocorr->level(level).stderror()));
        }
        impl_->binning = table;
    }
    for (const std::vector<double> &level : *impl_->binning)
        errors.push_back(level.at(i));
    return errors;
}

bool observable_result::is_scalar() const { return !impl_->declaration.vector_ && impl_->size <= 1; }

bool observable_result::is_histogram() const { return impl_->declaration.histogram_; }

const std::vector<std::uint64_t> &observable_result::histogram() const
{
    if (!impl_->declaration.histogram_)
        throw std::logic_error(quoted(impl_->declaration.name()) + " is not a histogram");
    return impl_->counts;
}

std::int64_t observable_result::histogram_min() const { return impl_->declaration.min_; }

const std::string &observable_result::sign() const { return impl_->declaration.sign_; }

observable_result transform(const vector_function &f, const std::vector<observable_result> &operands)
{
    if (operands.empty())
        throw std::invalid_argument("A function of results needs at least one observable_result");

    // The function of the parts of all operands: each operand first becomes
    // its own means, then `f` is applied to those.
    std::shared_ptr<observable_result::impl> derived(new observable_result::impl(observable("")));
    std::vector<vector_function> inner;
    std::vector<size_t> inner_size;
    for (const observable_result &operand : operands) {
        const observable_result::impl &source = *operand.impl_;
        source.need_bins();
        derived->parts.insert(derived->parts.end(), source.parts.begin(), source.parts.end());
        for (error_convergence component : operand.converged_vector())
            derived->convergence = std::max(derived->convergence, component);
        inner.push_back(source.function);
        inner_size.push_back(source.parts_size());
    }
    derived->function = [f, inner, inner_size](const std::vector<double> &means) {
        std::vector<double> arguments;
        std::vector<double>::const_iterator begin = means.begin();
        for (size_t i = 0; i != inner.size(); ++i) {
            const std::vector<double> part(begin, begin + inner_size[i]);
            const std::vector<double> value = inner[i] ? inner[i](part) : part;
            arguments.insert(arguments.end(), value.begin(), value.end());
            begin += inner_size[i];
        }
        return f(arguments);
    };
    derived->count = operands.front().count();
    derived->size = derived->evaluate().size();

    observable_result out;
    out.impl_ = derived;
    return out;
}

namespace {

observable_result binary(const observable_result &a, const observable_result &b, double (*op)(double, double))
{
    const size_t na = a.size(), nb = b.size();
    if (na != nb && na != 1 && nb != 1) {
        throw std::invalid_argument("Results with " + std::to_string(na) + " and " + std::to_string(nb)
                                    + " components cannot be combined");
    }
    observable_result out = transform([na, nb, op](const std::vector<double> &x) {
        std::vector<double> out(std::max(na, nb));
        for (size_t i = 0; i != out.size(); ++i)
            out[i] = op(x[na == 1 ? 0 : i], x[na + (nb == 1 ? 0 : i)]);
        return out;
    }, {a, b});
    if (!a.is_scalar() || !b.is_scalar())
        observables_access::set_vector(out);
    return out;
}

double add(double x, double y) { return x + y; }
double subtract(double x, double y) { return x - y; }
double multiply(double x, double y) { return x * y; }
double divide(double x, double y) { return x / y; }

}

observable_result apply(const std::function<double(double)> &f, const observable_result &a)
{
    observable_result out = transform([f](const std::vector<double> &x) {
        std::vector<double> out(x.size());
        std::transform(x.begin(), x.end(), out.begin(), f);
        return out;
    }, {a});
    if (!a.is_scalar())
        observables_access::set_vector(out);
    return out;
}

observable_result operator+(const observable_result &a, const observable_result &b) { return binary(a, b, add); }
observable_result operator-(const observable_result &a, const observable_result &b) { return binary(a, b, subtract); }
observable_result operator*(const observable_result &a, const observable_result &b) { return binary(a, b, multiply); }
observable_result operator/(const observable_result &a, const observable_result &b) { return binary(a, b, divide); }
observable_result operator+(const observable_result &a, double b) { return apply([b](double x) { return x + b; }, a); }
observable_result operator-(const observable_result &a, double b) { return apply([b](double x) { return x - b; }, a); }
observable_result operator*(const observable_result &a, double b) { return apply([b](double x) { return x * b; }, a); }
observable_result operator/(const observable_result &a, double b) { return apply([b](double x) { return x / b; }, a); }
observable_result operator+(double a, const observable_result &b) { return apply([a](double x) { return a + x; }, b); }
observable_result operator-(double a, const observable_result &b) { return apply([a](double x) { return a - x; }, b); }
observable_result operator*(double a, const observable_result &b) { return apply([a](double x) { return a * x; }, b); }
observable_result operator/(double a, const observable_result &b) { return apply([a](double x) { return a / x; }, b); }
observable_result operator-(const observable_result &a) { return apply([](double x) { return -x; }, a); }

// ACCUMULATOR

accumulator::accumulator(const observable &declaration) : impl_(new impl(declaration)) { }

void accumulator::add(const double *value, std::size_t size)
{
    impl &self = *impl_;
    if (self.declaration.histogram_) {
        const double offset = *value - self.declaration.min_;
        if (size != 1 || offset != std::floor(offset) || offset < 0 || offset >= self.counts.size()) {
            throw std::out_of_range("The histogram " + quoted(self.declaration.name())
                                    + " counts single integers within its range");
        }
        ++self.counts[size_t(offset)];
        return;
    }
    self.set_size(size);
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > sample(value, size);
    if (self.simple) {
        self.simple->add_sample(value);
    } else {
        *self.autocorr << make_adapter(sample);
        *self.batch << make_adapter(sample);
    }
}

void accumulator::add_vector(const double *value, std::size_t size)
{
    add(value, size);
    if (impl_->declaration.components_ != 1)
        impl_->declaration.vector_ = true;
}

std::uint64_t accumulator::count() const { return impl_->count(); }

const std::string &accumulator::name() const { return impl_->declaration.name(); }

observable_result accumulator::result() const
{
    if (!impl_->declaration.sign_.empty()) {
        throw std::logic_error(quoted(name()) + " is signed by " + quoted(impl_->declaration.sign_)
                               + "; its result comes from the results of its set");
    }
    observable_result out;
    out.impl_ = impl_->make_result();
    return out;
}

void accumulator::reset() { impl_->reset(); }

// SETS

observable_set::observable_set() : impl_(new impl) { }

observable_set::observable_set(const observable_set &other) : impl_(new impl) { *this = other; }

observable_set &observable_set::operator=(const observable_set &other)
{
    if (this == &other)
        return *this;
    impl_->observables.clear();
    for (const auto &named : other.impl_->observables) {
        accumulator copy(named.second);
        copy.impl_.reset(new accumulator::impl(*named.second.impl_));
        impl_->observables.insert(std::make_pair(named.first, copy));
    }
    return *this;
}

observable_set::~observable_set() { }

observable_set &observable_set::operator<<(const observable &declaration)
{
    if (declaration.name().empty())
        throw std::invalid_argument("An observable needs a name");
    if (!impl_->observables.insert(std::make_pair(declaration.name(), accumulator(declaration))).second)
        throw std::invalid_argument("The observable " + quoted(declaration.name()) + " is declared twice");
    return *this;
}

bool observable_set::has(const std::string &name) const { return impl_->observables.count(name) != 0; }

accumulator &observable_set::operator[](const std::string &name)
{
    const std::map<std::string, accumulator>::iterator found = impl_->observables.find(name);
    if (found == impl_->observables.end())
        throw std::out_of_range("There is no observable " + quoted(name));
    return found->second;
}

const accumulator &observable_set::operator[](const std::string &name) const
{
    return const_cast<observable_set &>(*this)[name];
}

std::vector<std::string> observable_set::names() const
{
    std::vector<std::string> names;
    for (const auto &named : impl_->observables)
        names.push_back(named.first);
    return names;
}

void observable_set::insert(const std::string &name, const accumulator &observable)
{
    if (name.empty())
        throw std::invalid_argument("An observable needs a name");
    impl_->observables.erase(name);
    impl_->observables.insert(std::make_pair(name, observable));
}

void observable_set::erase(const std::string &name)
{
    if (impl_->observables.erase(name) == 0)
        throw std::out_of_range("There is no observable " + quoted(name));
}

void observable_set::reset()
{
    for (auto &named : impl_->observables)
        named.second.reset();
}

result_set observable_set::results() const
{
    std::map<std::string, std::shared_ptr<observable_result::impl> > parts;
    for (const auto &named : impl_->observables) {
        parts[named.first] = named.second.impl_->make_result();
        parts[named.first]->declaration.name_ = named.first;
    }
    result_set out;
    out.impl_->link(parts);
    return out;
}

void observable_set::merge(const result_set &other)
{
    for (const auto &named : other.impl_->results) {
        const observable_result::impl &theirs = *named.second.impl_;
        if (!has(named.first))
            *this << theirs.declaration;
        accumulator::impl &mine = *(*this)[named.first].impl_;
        if (mine.declaration.histogram_ != theirs.declaration.histogram_
            || mine.declaration.with_bins_ != theirs.declaration.with_bins_
            || mine.declaration.sign_ != theirs.declaration.sign_) {
            throw std::invalid_argument(quoted(named.first) + " is declared differently in the two runs");
        }
        if (theirs.count == 0)
            continue;
        mine.declaration.vector_ = mine.declaration.vector_ || theirs.declaration.vector_;
        if (theirs.stored) {
            throw std::invalid_argument(
                quoted(named.first) + " was read from a result file. Runs are merged from their"
                " checkpoints or before they are written, not from stored results.");
        }
        if (theirs.declaration.histogram_) {
            if (mine.counts.size() != theirs.counts.size() || mine.declaration.min_ != theirs.declaration.min_)
                throw std::invalid_argument("The histogram " + quoted(named.first) + " differs in range");
            for (size_t i = 0; i != mine.counts.size(); ++i)
                mine.counts[i] += theirs.counts[i];
            continue;
        }
        mine.set_size(theirs.size);
        if (mine.simple) {
            *mine.simple << *theirs.simple;
        } else {
            *mine.autocorr << *theirs.autocorr;
            *mine.batch << *theirs.batch;
        }
    }
}

result_set::result_set() : impl_(new impl) { }
result_set::result_set(const result_set &other) : impl_(new impl(*other.impl_)) { }
result_set &result_set::operator=(const result_set &other) { *impl_ = *other.impl_; return *this; }
result_set::~result_set() { }

bool result_set::has(const std::string &name) const { return impl_->results.count(name) != 0; }

const observable_result &result_set::operator[](const std::string &name) const
{
    const std::map<std::string, observable_result>::const_iterator found = impl_->results.find(name);
    if (found == impl_->results.end())
        throw std::out_of_range("There is no result " + quoted(name));
    return found->second;
}

std::vector<std::string> result_set::names() const
{
    std::vector<std::string> names;
    for (const auto &named : impl_->results)
        names.push_back(named.first);
    return names;
}

void result_set::insert(const std::string &name, const observable_result &result)
{
    impl_->results[name] = result;
}

void result_set::erase(const std::string &name)
{
    if (impl_->results.erase(name) == 0)
        throw std::out_of_range("There is no result " + quoted(name));
}

result_set result_set::select(const std::vector<std::string> &names) const
{
    result_set out;
    for (const std::string &name : names)
        out.impl_->results[name] = (*this)[name];
    return out;
}

std::ostream &operator<<(std::ostream &out, const observable_result &result)
{
    if (result.count() == 0)
        return out << "no measurements";
    std::vector<double> mean, error, tau;
    std::vector<error_convergence> converged;
    try {
        mean = result.mean_vector();
        error = result.error_vector();
        tau = result.tau_vector();
        converged = result.converged_vector();
    } catch (const std::logic_error &) {
        const std::vector<std::uint64_t> &counts = result.histogram();
        out << "[";
        for (size_t i = 0; i != counts.size(); ++i)
            out << (i ? ", " : "") << counts[i];
        return out << "]";
    }
    if (mean.size() != 1)
        out << "[";
    for (size_t i = 0; i != mean.size(); ++i) {
        out << (i ? ", " : "") << mean[i] << " +/- " << error[i];
        if (tau[i] == tau[i])
            out << " (tau = " << tau[i] << ")";
        if (converged[i] == MAYBE_CONVERGED)
            out << " WARNING: check error convergence";
        if (converged[i] == NOT_CONVERGED)
            out << " WARNING: ERRORS NOT CONVERGED!!!";
    }
    if (mean.size() != 1)
        out << "]";
    return out;
}

std::ostream &operator<<(std::ostream &out, const result_set &results)
{
    for (const std::string &name : results.names())
        out << name << ": " << results[name] << "\n";
    return out;
}

namespace {

/**
 * Groups the small reductions of a whole result set into one per type until
 * `commit`. The interface allows a reducer to defer; this one does so for any
 * reducer it wraps, without holding a copy of the large ones: those are
 * passed on at once, as their cost is the data and not the latency.
 */
class grouping_reducer : public reducer
{
public:
    explicit grouping_reducer(const reducer &inner) : inner_(inner) { }

    reducer_setup get_setup() const override { return inner_.get_setup(); }
    int64_t get_max(int64_t value) const override { return inner_.get_max(value); }
    void reduce(view<double> data) const override { defer(data, doubles_); }
    void reduce(view<int32_t> data) const override { defer(data, ints_); }
    void reduce(view<int64_t> data) const override { defer(data, longs_); }

    void commit() const override
    {
        flush(doubles_);
        flush(ints_);
        flush(longs_);
        inner_.commit();
    }

private:
    static const size_t small = 4096;

    template <typename T>
    void defer(view<T> data, std::vector<view<T> > &pending) const
    {
        if (data.size() > small)
            inner_.reduce(data);
        else if (data.size() != 0)
            pending.push_back(data);
    }

    template <typename T>
    void flush(std::vector<view<T> > &pending) const
    {
        size_t total = 0;
        for (const view<T> &data : pending)
            total += data.size();
        std::vector<T> packed;
        packed.reserve(total);
        for (const view<T> &data : pending)
            packed.insert(packed.end(), data.data(), data.data() + data.size());
        if (total != 0)
            inner_.reduce(view<T>(packed.data(), packed.size()));
        size_t offset = 0;
        for (view<T> &data : pending) {
            std::copy(packed.begin() + offset, packed.begin() + offset + data.size(), data.data());
            offset += data.size();
        }
        pending.clear();
    }

    const reducer &inner_;
    mutable std::vector<view<double> > doubles_;
    mutable std::vector<view<int32_t> > ints_;
    mutable std::vector<view<int64_t> > longs_;
};

}

namespace {

// The result itself if nothing but `held` refers to it and `held` belongs to
// one set only, else a copy. For an observable of a gigabyte, a copy costs a
// gigabyte more on every rank.
template <typename R>
std::unique_ptr<R> take_or_copy(std::shared_ptr<const R> &held, bool sole)
{
    if (sole && held.use_count() == 1) {
        // Results are created mutable and only held as const.
        std::unique_ptr<R> out(new R(std::move(const_cast<R &>(*held))));
        held.reset();
        return out;
    }
    return std::unique_ptr<R>(new R(*held));
}

}

void result_set::reduce(const reducer &outer)
{
    // All observables first post their reductions, then one commit sends the
    // small ones together, then each collects its result.
    const grouping_reducer r(outer);
    const reducer_setup setup = r.get_setup();
    struct pending
    {
        std::shared_ptr<observable_result::impl> reduced;
        std::vector<int64_t> counts;
        std::unique_ptr<autocorr_result<double> > autocorr;
        std::unique_ptr<batch_result<double> > batch;
        std::unique_ptr<var_result<double> > simple;
    };
    std::map<std::string, std::shared_ptr<observable_result::impl> > parts;
    std::vector<pending> all;
    all.reserve(impl_->results.size());
    // If a result of this set was changed for the reduction and the reduction
    // fails, the set is left empty rather than with a part of its data.
    bool changed = false;
    struct empty_on_failure
    {
        result_set::impl &set;
        const bool &changed;
        bool done;
        ~empty_on_failure() { if (!done && changed) set.results.clear(); }
    } guard{*impl_, changed, false};
    for (const auto &named : impl_->results) {
        // Results are created mutable and only held as const.
        observable_result::impl &mine = const_cast<observable_result::impl &>(*named.second.impl_);
        if (mine.stored)
            throw std::invalid_argument(quoted(named.first) + " was read from a result file and cannot be reduced");
        // A result no other set shares is replaced by its reduction below, so
        // its data can be reduced in place. Its own references to bins go
        // first, those of a signed observable to the bins of its sign too;
        // data that anything else still refers to is copied.
        const bool sole = named.second.impl_.use_count() == 1;
        if (sole) {
            changed = true;
            mine.parts.clear();
            mine.evaluated.reset();
            mine.binning.reset();
        }
        all.emplace_back();
        pending &p = all.back();
        p.reduced.reset(new observable_result::impl(mine.declaration));
        parts[named.first] = p.reduced;

        if (mine.declaration.histogram_) {
            p.counts.assign(mine.counts.begin(), mine.counts.end());
            r.reduce(view<int64_t>(p.counts.data(), p.counts.size()));
            continue;
        }

        // An instance without samples does not know the number of components.
        p.reduced->size = r.get_max(mine.size);
        if (p.reduced->size == 0)
            continue;
        if (mine.size != 0 && mine.size != p.reduced->size)
            throw std::invalid_argument(quoted(named.first) + " differs in size between instances");
        if (mine.declaration.with_bins_) {
            if (mine.autocorr)
                p.autocorr = take_or_copy(mine.autocorr, sole);
            else
                p.autocorr.reset(new autocorr_result<double>(autocorr_acc<double>(p.reduced->size).finalize()));
            if (mine.batch)
                p.batch = take_or_copy(mine.batch, sole);
            else
                p.batch.reset(new batch_result<double>(
                    batch_acc<double>(p.reduced->size, mine.declaration.num_bins_).finalize()));
            p.autocorr->reduce(r, true, false);
            p.batch->reduce(r, true, false);
        } else {
            if (mine.simple)
                p.simple = take_or_copy(mine.simple, sole);
            else
                p.simple.reset(new var_result<double>(var_acc<double>(p.reduced->size).finalize()));
            p.simple->reduce(r, true, false);
        }
    }
    r.commit();
    for (pending &p : all) {
        observable_result::impl &reduced = *p.reduced;
        if (reduced.declaration.histogram_) {
            reduced.counts.assign(p.counts.begin(), p.counts.end());
            reduced.size = p.counts.size();
            reduced.count = std::accumulate(p.counts.begin(), p.counts.end(), uint64_t(0));
        } else if (p.autocorr) {
            p.autocorr->reduce(r, false, true);
            p.batch->reduce(r, false, true);
            if (setup.have_result) {
                reduced.count = p.autocorr->count();
                reduced.autocorr.reset(p.autocorr.release());
                reduced.batch.reset(p.batch.release());
                reduced.parts.push_back(reduced.batch);
            }
        } else if (p.simple) {
            p.simple->reduce(r, false, true);
            if (setup.have_result) {
                reduced.count = p.simple->count();
                reduced.simple.reset(p.simple.release());
            }
        }
    }
    if (setup.have_result)
        impl_->link(parts);
    else
        impl_->results.clear();
    guard.done = true;
}

}} /* namespace alps::alea */
