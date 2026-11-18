/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

/**
 * Named observables of a Monte Carlo simulation.
 *
 * This is what a simulation and a program that evaluates one work with:
 *
 *     alps::alea::observable_set measurements;
 *     measurements << observable("Energy")
 *                  << observable("Sign")
 *                  << observable("Density").signed_by("Sign")
 *                  << observable("Counter").without_bins();
 *
 *     measurements["Energy"] << energy;               // once per sweep
 *     measurements["Sign"] << sign;
 *     measurements["Density"] << sign * density;      // a vector
 *
 *     alps::alea::result_set results = measurements.results();
 *     double e = results["Energy"].mean(), de = results["Energy"].error();
 *     alps::alea::result c = results["Energy"] * results["Energy"];
 *
 * An observable keeps a binning analysis, for its error and autocorrelation
 * time, and a fixed number of bins of equal length in time order, for
 * functions of several observables. Such functions are evaluated by jackknife
 * over the bins and need bins on every operand that line up; there is no
 * fallback that treats operands as independent.
 *
 * The accumulators of `<alps/alea.hpp>` do the work. This header does not
 * expose them, or Eigen.
 */
#pragma once

#include <alps/alea/error_convergence.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace alps { namespace alea {

struct reducer;

class observable_result;
class result_set;
class observable_set;
struct observables_access;

/**
 * Declaration of an observable, for `observable_set << observable(name)`.
 *
 * The number of components is that of the first sample.
 */
class observable
{
public:
    explicit observable(const std::string &name);

    /**
     * Keep mean and naive error only, for samples known to be uncorrelated.
     * Such an observable cannot enter a function of observables.
     */
    observable &without_bins();

    /** Number of bins kept for jackknife; even, 128 unless set */
    observable &bins(std::size_t num_bins);

    /**
     * The samples are the quantity times the sign named here; the result is
     * the ratio of the two averages, evaluated by jackknife. The sign must be
     * an observable of the same set, measured in the same sweeps.
     */
    observable &signed_by(const std::string &sign);

    /** Names of the components */
    observable &labels(const std::vector<std::string> &labels);

    const std::string &name() const { return name_; }

private:
    friend struct observables_access;
    friend class accumulator;
    friend class observable_set;
    friend class observable_result;
    friend class result_set;
    std::string name_, sign_;
    std::vector<std::string> labels_;
    std::size_t num_bins_;
    bool with_bins_, histogram_;
    std::int64_t min_, max_;
    friend observable histogram(const std::string &, std::int64_t, std::int64_t);
};

/** Declaration of a histogram of the integers `min` to `max`, inclusive */
observable histogram(const std::string &name, std::int64_t min, std::int64_t max);

/**
 * One observable of a set, while it is being measured.
 */
class accumulator
{
public:
    /** Add a sample. A histogram counts the integer value given. */
    accumulator &operator<<(double value) { add(&value, 1); return *this; }

    accumulator &operator<<(const std::vector<double> &value) { add(value.data(), value.size()); return *this; }

    void add(const double *value, std::size_t size);

    /** Number of samples */
    std::uint64_t count() const;

    /** Discard all samples */
    void reset();

    struct impl;    // not part of the interface

private:
    friend struct observables_access;
    friend class observable_set;
    explicit accumulator(const observable &declaration);
    std::shared_ptr<impl> impl_;
};

/**
 * The observables of a simulation, by name.
 */
class observable_set
{
public:
    observable_set();
    observable_set(const observable_set &other);
    observable_set &operator=(const observable_set &other);
    ~observable_set();

    /** Declare an observable. Declaring a name twice is an error. */
    observable_set &operator<<(const observable &declaration);

    bool has(const std::string &name) const;

    /** The observable of that name; an unknown name is an error */
    accumulator &operator[](const std::string &name);
    const accumulator &operator[](const std::string &name) const;

    /** Names in alphabetical order */
    std::vector<std::string> names() const;

    /** Discard all samples and keep the declarations, as after thermalization */
    void reset();

    /** What has been measured so far; the set can go on measuring */
    result_set results() const;

    /**
     * Append the measurements of an independent run to these: counts add up,
     * bins are appended. Observables this set does not have are added.
     */
    void merge(const result_set &other);

private:
    friend struct observables_access;
    struct impl;
    std::unique_ptr<impl> impl_;
};

/**
 * An observable after measuring, or a function of such.
 *
 * `mean()`, `error()` and `tau()` are for observables with one component and
 * refuse others; the `_vector` forms work for any.
 */
class observable_result
{
public:
    observable_result();
    observable_result(const observable_result &other);
    observable_result &operator=(const observable_result &other);
    ~observable_result();

    /** Number of components; 0 before the first sample */
    std::size_t size() const;

    /** Number of samples */
    std::uint64_t count() const;

    double mean() const;
    double error() const;
    double tau() const;
    error_convergence converged() const;

    std::vector<double> mean_vector() const;
    std::vector<double> error_vector() const;

    /** Variance of a single sample */
    std::vector<double> variance_vector() const;

    /** Integrated autocorrelation time; NaN where it is not known */
    std::vector<double> tau_vector() const;

    std::vector<error_convergence> converged_vector() const;

    const std::vector<std::string> &labels() const;

    /** Can this result enter a function of results? */
    bool has_bins() const;

    /** Number of samples in a full bin */
    std::uint64_t bin_size() const;

    /** Means of the full bins in time order, one vector per bin */
    std::vector<std::vector<double> > bin_means() const;

    /**
     * The binning analysis: the error of the mean of component `i` estimated
     * from bins of 1, 2, 4, ... samples. Empty without a binning analysis.
     */
    std::vector<double> binning_errors(std::size_t i=0) const;

    /** Counts of a histogram, from its minimum to its maximum */
    const std::vector<std::uint64_t> &histogram() const;
    std::int64_t histogram_min() const;

    /** Name of the sign of a signed observable, or empty */
    const std::string &sign() const;

    struct impl;    // not part of the interface

private:
    friend class observable_set;
    friend class result_set;
    friend class accumulator;
    friend observable_result transform(const std::function<std::vector<double>(const std::vector<double> &)> &,
                            const std::vector<observable_result> &);
    friend struct observables_access;
    std::shared_ptr<const impl> impl_;
};

/**
 * A function of results by jackknife: `f` maps the means of all operands, one
 * after the other, to the means of the result. Every operand needs bins, and
 * the bins must hold the same numbers of samples; otherwise this throws
 * `std::invalid_argument`.
 */
observable_result transform(const std::function<std::vector<double>(const std::vector<double> &)> &f,
                 const std::vector<observable_result> &operands);

/** Component by component; a result with one component combines with any */
observable_result operator+(const observable_result &a, const observable_result &b);
observable_result operator-(const observable_result &a, const observable_result &b);
observable_result operator*(const observable_result &a, const observable_result &b);
observable_result operator/(const observable_result &a, const observable_result &b);
observable_result operator+(const observable_result &a, double b);
observable_result operator-(const observable_result &a, double b);
observable_result operator*(const observable_result &a, double b);
observable_result operator/(const observable_result &a, double b);
observable_result operator+(double a, const observable_result &b);
observable_result operator-(double a, const observable_result &b);
observable_result operator*(double a, const observable_result &b);
observable_result operator/(double a, const observable_result &b);
observable_result operator-(const observable_result &a);

/** `f` applied to every component */
observable_result apply(const std::function<double(double)> &f, const observable_result &a);

/**
 * The results of a simulation, by name.
 */
class result_set
{
public:
    result_set();
    result_set(const result_set &other);
    result_set &operator=(const result_set &other);
    ~result_set();

    bool has(const std::string &name) const;

    /** The result of that name; an unknown name is an error */
    const observable_result &operator[](const std::string &name) const;

    std::vector<std::string> names() const;

    /**
     * Collect the results of all instances of a parallel run, as `merge` does
     * for one run after the other. All instances must have declared the same
     * observables. Instances without the result are left with an empty set.
     */
    void reduce(const reducer &r);

private:
    friend struct observables_access;
    friend class observable_set;
    struct impl;
    std::unique_ptr<impl> impl_;
};

}} /* namespace alps::alea */
