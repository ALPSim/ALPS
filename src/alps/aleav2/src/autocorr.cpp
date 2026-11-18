/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include <alps/alea/autocorr.hpp>
#include <alps/alea/serialize.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/state.hpp>
#include <alps/alea/internal/format.hpp>

#include <cmath>
#include <limits>

namespace alps { namespace alea {

template <typename T>
autocorr_acc<T>::autocorr_acc(size_t size, uint64_t batch_size, size_t granularity)
    : size_(size)
    , batch_size_(batch_size)
    , count_(0)
    , nextlevel_(batch_size)
    , granularity_(granularity)
    , level_()
{
    level_.push_back(var_acc<T>(size, batch_size));
}

template <typename T>
void autocorr_acc<T>::reset()
{
    count_ = 0;
    nextlevel_ = batch_size_;
    level_.clear();
    level_.push_back(var_acc<T>(size_, batch_size_));
}

template <typename T>
void autocorr_acc<T>::set_size(size_t size)
{
    size_ = size;
    reset();
}

template <typename T>
void autocorr_acc<T>::set_batch_size(uint64_t batch_size)
{
    // TODO: handle the case where we just discard levels more gracefully
    batch_size_ = batch_size;
    reset();
}

template <typename T>
void autocorr_acc<T>::set_granularity(size_t granularity)
{
    // TODO: handle the case where we just discard levels more gracefully
    granularity_ = granularity;
    reset();
}

template <typename T>
void autocorr_acc<T>::add_level()
{
    // add a new level on top and push back the nextlevel
    nextlevel_ *= granularity_;
    level_.push_back(var_acc<T>(size(), nextlevel_));
}

template <typename T>
void autocorr_acc<T>::add(const computed<T> &source, uint64_t count)
{
    assert(count_ < nextlevel_);
    internal::check_valid(*this);

    // if we require next level, then do it!
    count_ += count;
    if(count_ >= nextlevel_)
        add_level();

    // now add current element at the bottom and watch it propagate
    level_[0].add(source, count, level_.data() + 1);
}

template <typename T>
autocorr_acc<T> &autocorr_acc<T>::operator<<(const autocorr_result<T> &other)
{
    internal::check_valid(*this);
    internal::check_valid(other);
    if (size() != other.size())
        throw size_mismatch();

    // Levels that both have are merged batch by batch. Where this accumulator
    // is the shorter one, its own samples reach the added levels as they do in
    // a single run: when a batch below fills up, or when the result is taken.
    const size_t other_nlevel = other.nlevel();
    while (nlevel() < other_nlevel)
        add_level();
    for (size_t i = 0; i != other_nlevel; ++i)
        level_[i] << other.level(i);

    // As in a single run, the highest level must be too large to fill up.
    count_ += other.count();
    while (count_ >= nextlevel_)
        add_level();

    // Levels above those of the result know nothing of its samples yet. They
    // amount to less than one batch there, and enter as such.
    if (nlevel() > other_nlevel) {
        const column<T> sum = other.mean() * double(other.count());
        level_[other_nlevel].add(make_adapter(sum), other.count(),
                                 level_.data() + other_nlevel + 1);
    }
    return *this;
}

template <typename T>
autocorr_result<T> autocorr_acc<T>::result() const
{
    internal::check_valid(*this);
    autocorr_result<T> result;
    autocorr_acc<T>(*this).finalize_to(result);
    return result;
}

template <typename T>
autocorr_result<T> autocorr_acc<T>::finalize()
{
    autocorr_result<T> result;
    finalize_to(result);
    return result;
}

template <typename T>
void autocorr_acc<T>::finalize_to(autocorr_result<T> &result)
{
    internal::check_valid(*this);
    result.level_.resize(level_.size());

    // Finalize each level.
    // NOTE: it is imperative to do this in bottom-up order, since it collects
    //       the left-over data in the current batch from the low-lying levels
    //       and propagates them upwards.
    for (size_t i = 0; i != level_.size() - 1; ++i)
        level_[i].finalize_to(result.level_[i], level_.data() + i + 1);

    level_[nlevel() - 1].finalize_to(result.level_[nlevel() - 1], nullptr);
    level_.clear();     // signal invalidity
}

template <typename T>
void autocorr_acc<T>::save_state(serializer &s, const std::string &key) const
{
    internal::check_valid(*this);
    internal::serializer_sentry group(s, key);
    internal::write_state_header(s, internal::STATE_AUTOCORR);
    serialize(s, "@size", static_cast<uint64_t>(size_));
    serialize(s, "@batch_size", static_cast<uint64_t>(batch_size_));
    serialize(s, "@granularity", static_cast<uint64_t>(granularity_));
    serialize(s, "@nlevel", static_cast<uint64_t>(nlevel()));
    serialize(s, "count", static_cast<uint64_t>(count_));
    serialize(s, "nextlevel", static_cast<uint64_t>(nextlevel_));

    s.enter("level");
    for (size_t i = 0; i != nlevel(); ++i)
        level_[i].save_state(s, std::to_string(i));
    s.exit();
}

template <typename T>
void autocorr_acc<T>::load_state(deserializer &s, const std::string &key)
{
    internal::deserializer_sentry group(s, key);
    internal::read_state_header(s, internal::STATE_AUTOCORR);
    size_ = internal::read_uint(s, "@size");
    batch_size_ = internal::read_uint(s, "@batch_size");
    granularity_ = internal::read_uint(s, "@granularity");
    const size_t new_nlevel = internal::read_uint(s, "@nlevel");
    count_ = internal::read_uint(s, "count");
    nextlevel_ = internal::read_uint(s, "nextlevel");

    level_.assign(new_nlevel, level_acc_type(size_));
    s.enter("level");
    for (size_t i = 0; i != new_nlevel; ++i)
        level_[i].load_state(s, std::to_string(i));
    s.exit();
}

template class autocorr_acc<double>;
template class autocorr_acc<std::complex<double> >;

template <typename T>
bool operator==(const autocorr_result<T> &r1, const autocorr_result<T> &r2)
{
    if(r1.nlevel() != r2.nlevel()) return false;
    for(size_t i = 0; i < r1.nlevel(); ++i) {
        if(r1.level(i) != r2.level(i))
            return false;
    }
    return true;
}

template bool operator==(const autocorr_result<double> &r1,
                         const autocorr_result<double> &r2);
template bool operator==(const autocorr_result<std::complex<double>> &r1,
                         const autocorr_result<std::complex<double>> &r2);

template <typename T>
uint64_t autocorr_result<T>::batch_size(size_t i) const
{
    return level_[i].batch_size();
}

template <typename T>
size_t autocorr_result<T>::find_level(uint64_t min_batches) const
{
    for (size_t i = nlevel(); i != 0; --i) {
        if (level(i-1).count() / level(i-1).batch_size() >= min_batches)
            return i - 1;
    }
    return 0;
}

template <typename T>
double autocorr_result<T>::count2(uint64_t min_batches) const
{
    return level_[find_level(min_batches)].count2();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::var(uint64_t min_batches) const
{
    // The question is whether to return the ensemble variance or the batch
    // variance here.  We opt for the batch variance, because it makes
    // everything consistent with each other.
    return level_[find_level(min_batches)].var();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::stderror(uint64_t min_batches) const
{
    return level_[find_level(min_batches)].stderror();
}

template <typename T>
double autocorr_result<T>::observations(uint64_t min_batches) const
{
    return level_[find_level(min_batches)].observations();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::tau(uint64_t min_batches) const
{
    const size_t lvl = find_level(min_batches);
    if (lvl == 0)
        return column<var_type>::Constant(size(), std::numeric_limits<double>::quiet_NaN());

    const column<var_type> &var0 = level_[0].var();
    const column<var_type> &varn = level_[lvl].var();
    const double n = double(level_[lvl].batch_size()) / level_[0].batch_size();

    return (0.5 * n * varn.array() / var0.array() - 0.5).matrix();
}

template <typename T>
std::vector<error_convergence> autocorr_result<T>::converged(uint64_t min_batches) const
{
    const size_t range = 4;
    const size_t lvl = find_level(min_batches);
    if (lvl + 1 < range)
        return std::vector<error_convergence>(size(), MAYBE_CONVERGED);

    std::vector<error_convergence> conv(size(), CONVERGED);
    const column<var_type> err = level_[lvl].stderror();
    for (size_t i = lvl + 1 - range; i != lvl; ++i) {
        const column<var_type> this_err = level_[i].stderror();
        for (size_t k = 0; k != size(); ++k) {
            const double ratio = std::abs(this_err[k]), full = std::abs(err[k]);
            if (ratio >= full)
                conv[k] = CONVERGED;
            else if (ratio < 0.824 * full)
                conv[k] = NOT_CONVERGED;
            else if (ratio < 0.9 * full && conv[k] != NOT_CONVERGED)
                conv[k] = MAYBE_CONVERGED;
        }
    }
    return conv;
}

template <typename T>
void autocorr_result<T>::reduce(const reducer &r, bool pre_commit, bool post_commit)
{
    internal::check_valid(*this);

    if (pre_commit) {
        // Bring all instances to the same number of levels. To the levels a
        // shorter run lacks, its samples are less than one batch.
        const size_t needs_levels = r.get_max(nlevel());
        if (nlevel() < needs_levels) {
            var_data<T> one_batch(size());
            if (count() != 0) {
                one_batch.data() = mean();
                one_batch.data2().fill(0);
                one_batch.count() = count();
                one_batch.count2() = double(count()) * count();
            }
            const level_result_type all(one_batch);
            level_.resize(needs_levels, all);
        }

        for (size_t i = 0; i != nlevel(); ++i)
            level_[i].reduce(r, true, false);
    }
    if (pre_commit && post_commit) {
        // perform commit
        r.commit();
    }
    if (post_commit) {
        // cleanups
        reducer_setup setup = r.get_setup();
        for (size_t i = 0; i != nlevel(); ++i)
            level_[i].reduce(r, false, true);
        if (!setup.have_result)
            level_.clear();         // invalidate
    }
}

template class autocorr_result<double>;
template class autocorr_result<std::complex<double> >;


template <typename T>
void serialize(serializer &s, const std::string &key, const autocorr_result<T> &self)
{
    internal::check_valid(self);
    internal::serializer_sentry group(s, key);

    // Write size and nlevel as 64-bit integers for consistency
    serialize(s, "@size", static_cast<uint64_t>(self.size()));
    serialize(s, "@nlevel", static_cast<uint64_t>(self.nlevel()));

    s.enter("level");
    for (size_t i = 0; i != self.nlevel(); ++i)
        serialize(s, std::to_string(i), self.level_[i]);
    s.exit();

    s.enter("mean");
    serialize(s, "value", self.mean());
    serialize(s, "error", self.stderror());
    s.exit();
}

template <typename T>
void deserialize(deserializer &s, const std::string &key, autocorr_result<T> &self)
{
    typedef typename autocorr_result<T>::var_type var_type;
    internal::deserializer_sentry group(s, key);

    size_t scalar_size = 1;
    s.read("@size", ndview<uint64_t>(nullptr, &scalar_size, 0)); // discard

    // first deserialize the fundamentals and make sure that the target fits
    uint64_t new_nlevel;
    deserialize(s, "@nlevel", new_nlevel);
    self.level_.resize(new_nlevel);

    s.enter("level");
    for (size_t i = 0; i != self.nlevel(); ++i)
        deserialize(s, std::to_string(i), self.level_[i]);
    s.exit();

    scalar_size = self.size();
    s.enter("mean");
    s.read("value", ndview<T>(nullptr, &scalar_size, 1)); // discard
    s.read("error", ndview<var_type>(nullptr, &scalar_size, 1)); // discard
    s.exit();
}

template void serialize(serializer &, const std::string &key, const autocorr_result<double> &);
template void serialize(serializer &, const std::string &key, const autocorr_result<std::complex<double>> &);

template void deserialize(deserializer &, const std::string &key, autocorr_result<double> &);
template void deserialize(deserializer &, const std::string &key, autocorr_result<std::complex<double> > &);

template <typename T>
std::ostream &operator<<(std::ostream &str, const autocorr_result<T> &self)
{
    internal::check_valid(self);
    internal::format_sentry sentry(str);
    verbosity verb = internal::get_format(str, PRINT_TERSE);

    if (verb == PRINT_VERBOSE)
        str << "<X> = ";
    str << self.mean() << " +- " << self.stderror();

    if (verb == PRINT_VERBOSE) {
        str << "\nLevels:" << PRINT_TERSE;
        for (const var_result<T> &curr : self.level_)
            str << "\n  " << curr;
    }
    return str;
}

template std::ostream &operator<<(std::ostream &, const autocorr_result<double> &);
template std::ostream &operator<<(std::ostream &, const autocorr_result<std::complex<double>> &);

}}

