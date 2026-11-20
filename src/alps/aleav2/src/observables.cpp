/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include "observables_impl.hpp"

namespace alps { namespace alea {


// DECLARATION

observable::observable(const std::string &name)
    : name_(name)
    , num_bins_(128)
    , with_bins_(true)
    , histogram_(false)
    , min_(0)
    , max_(0)
{ }

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
    return transform([na, nb, op](const std::vector<double> &x) {
        std::vector<double> out(std::max(na, nb));
        for (size_t i = 0; i != out.size(); ++i)
            out[i] = op(x[na == 1 ? 0 : i], x[na + (nb == 1 ? 0 : i)]);
        return out;
    }, {a, b});
}

double add(double x, double y) { return x + y; }
double subtract(double x, double y) { return x - y; }
double multiply(double x, double y) { return x * y; }
double divide(double x, double y) { return x / y; }

}

observable_result apply(const std::function<double(double)> &f, const observable_result &a)
{
    return transform([f](const std::vector<double> &x) {
        std::vector<double> out(x.size());
        std::transform(x.begin(), x.end(), out.begin(), f);
        return out;
    }, {a});
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
        *self.simple << make_adapter(sample);
    } else {
        *self.autocorr << make_adapter(sample);
        *self.batch << make_adapter(sample);
    }
}

std::uint64_t accumulator::count() const { return impl_->count(); }

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

void observable_set::reset()
{
    for (auto &named : impl_->observables)
        named.second.reset();
}

result_set observable_set::results() const
{
    std::map<std::string, std::shared_ptr<observable_result::impl> > parts;
    for (const auto &named : impl_->observables)
        parts[named.first] = named.second.impl_->make_result();
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
        throw std::out_of_range("There is no observable_result " + quoted(name));
    return found->second;
}

std::vector<std::string> result_set::names() const
{
    std::vector<std::string> names;
    for (const auto &named : impl_->results)
        names.push_back(named.first);
    return names;
}

void result_set::reduce(const reducer &r)
{
    const reducer_setup setup = r.get_setup();
    std::map<std::string, std::shared_ptr<observable_result::impl> > parts;
    for (const auto &named : impl_->results) {
        const observable_result::impl &mine = *named.second.impl_;
        if (mine.stored)
            throw std::invalid_argument(quoted(named.first) + " was read from a result file and cannot be reduced");
        std::shared_ptr<observable_result::impl> reduced(new observable_result::impl(mine.declaration));
        parts[named.first] = reduced;

        if (mine.declaration.histogram_) {
            std::vector<int64_t> counts(mine.counts.begin(), mine.counts.end());
            r.reduce(view<int64_t>(counts.data(), counts.size()));
            r.commit();
            reduced->counts.assign(counts.begin(), counts.end());
            reduced->size = counts.size();
            reduced->count = std::accumulate(counts.begin(), counts.end(), uint64_t(0));
            continue;
        }

        // An instance without samples does not know the number of components.
        reduced->size = r.get_max(mine.size);
        if (reduced->size == 0)
            continue;
        if (mine.size != 0 && mine.size != reduced->size)
            throw std::invalid_argument(quoted(named.first) + " differs in size between instances");
        if (mine.declaration.with_bins_) {
            autocorr_result<double> autocorr = mine.autocorr
                ? *mine.autocorr : autocorr_acc<double>(reduced->size).finalize();
            batch_result<double> batch = mine.batch
                ? *mine.batch : batch_acc<double>(reduced->size, mine.declaration.num_bins_).finalize();
            autocorr.reduce(r);
            batch.reduce(r);
            if (setup.have_result) {
                reduced->count = autocorr.count();
                reduced->autocorr.reset(new autocorr_result<double>(autocorr));
                reduced->batch.reset(new batch_result<double>(batch));
                reduced->parts.push_back(reduced->batch);
            }
        } else {
            var_result<double> simple = mine.simple ? *mine.simple : var_acc<double>(reduced->size).finalize();
            simple.reduce(r);
            if (setup.have_result) {
                reduced->count = simple.count();
                reduced->simple.reset(new var_result<double>(simple));
            }
        }
    }
    if (setup.have_result)
        impl_->link(parts);
    else
        impl_->results.clear();
}

}} /* namespace alps::alea */
