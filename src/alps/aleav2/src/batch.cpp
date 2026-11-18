/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include <alps/alea/batch.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/serialize.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/format.hpp>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace alps { namespace alea {

namespace {

/**
 * The batches of a series as a list: full batches in order, and the samples
 * left over at the end pooled into one partial batch.
 */
template <typename T>
class batch_list
{
public:
    batch_list(size_t size, uint64_t batch_size)
        : batch_size_(batch_size)
        , partial_(column<T>::Zero(size))
        , partial_count_(0)
    { }

    /** Read the non-empty batches in the given range of columns */
    batch_list(const batch_data<T> &data, size_t begin, size_t end)
        : batch_list(data.size(), data.batch_size())
    {
        if (batch_size_ == 0) {
            // Not recorded: the smallest batch is a full one.
            for (size_t i = begin; i != end; ++i) {
                if (data.count()(i) != 0 && (batch_size_ == 0 || data.count()(i) < batch_size_))
                    batch_size_ = data.count()(i);
            }
            if (batch_size_ == 0)
                batch_size_ = 1;
        }
        for (size_t i = begin; i != end; ++i)
            add(data.batch().col(i), data.count()(i));
    }

    explicit batch_list(const batch_data<T> &data)
        : batch_list(data, 0, data.num_batches())
    { }

    uint64_t batch_size() const { return batch_size_; }

    size_t num_full() const { return full_.size(); }

    void add(const column<T> &sum, uint64_t count)
    {
        if (count >= batch_size_) {
            full_.push_back(sum);
            count_.push_back(count);
        } else if (count != 0) {
            partial_ += sum;
            partial_count_ += count;
            promote();
        }
    }

    void append(const batch_list &other)
    {
        for (size_t i = 0; i != other.full_.size(); ++i)
            add(other.full_[i], other.count_[i]);
        add(other.partial_, other.partial_count_);
    }

    /** Add neighbouring batches and double the batch size */
    void coarsen()
    {
        const size_t pairs = full_.size() / 2;
        for (size_t i = 0; i != pairs; ++i) {
            full_[i] = full_[2 * i] + full_[2 * i + 1];
            count_[i] = count_[2 * i] + count_[2 * i + 1];
        }
        if (full_.size() % 2 != 0) {
            partial_ += full_.back();
            partial_count_ += count_.back();
        }
        full_.resize(pairs);
        count_.resize(pairs);
        batch_size_ *= 2;
        promote();
    }

    void coarsen_to(uint64_t batch_size)
    {
        while (batch_size_ < batch_size)
            coarsen();
        if (batch_size_ != batch_size)
            throw std::invalid_argument("Batch sizes must differ by a power of two");
    }

    /** Coarsen until `store()` has room in the given number of batches */
    void fit(size_t num_batches)
    {
        while (full_.size() > num_batches
               || (full_.size() == num_batches && partial_count_ != 0))
            coarsen();
    }

    /** Write to columns starting at `begin`; returns the number of full batches */
    size_t store(batch_data<T> &data, size_t begin=0) const
    {
        for (size_t i = 0; i != full_.size(); ++i) {
            data.batch().col(begin + i) = full_[i];
            data.count()(begin + i) = count_[i];
        }
        if (partial_count_ != 0) {
            data.batch().col(begin + full_.size()) = partial_;
            data.count()(begin + full_.size()) = partial_count_;
        }
        data.set_batch_size(batch_size_);
        return full_.size();
    }

private:
    void promote()
    {
        if (partial_count_ < batch_size_)
            return;
        full_.push_back(partial_);
        count_.push_back(partial_count_);
        partial_.fill(0);
        partial_count_ = 0;
    }

    uint64_t batch_size_;
    std::vector<column<T> > full_;
    std::vector<uint64_t> count_;
    column<T> partial_;
    uint64_t partial_count_;
};

}

template <typename T>
batch_data<T>::batch_data(size_t size, size_t num_batches)
    : batch_(size, num_batches)
    , count_(num_batches)
    , batch_size_(0)
{
    reset();
}

template <typename T>
void batch_data<T>::reset()
{
    batch_.fill(0);
    count_.fill(0);
}

template class batch_data<double>;
template class batch_data<std::complex<double> >;


template <typename T>
batch_acc<T>::batch_acc(size_t size, size_t num_batches, uint64_t base_size)
    : size_(size)
    , num_batches_(num_batches)
    , base_size_(base_size)
    , store_(new batch_data<T>(size, num_batches))
    , cursor_(0)
{
    if (num_batches % 2 != 0) {
        throw std::runtime_error("Number of batches must be even to allow "
                                 "for rebatching.");
    }
    store_->set_batch_size(base_size_);
}

template <typename T>
batch_acc<T>::batch_acc(const batch_acc &other)
    : size_(other.size_)
    , num_batches_(other.num_batches_)
    , base_size_(other.base_size_)
    , store_(other.store_ ? new batch_data<T>(*other.store_) : nullptr)
    , cursor_(other.cursor_)
{ }

template <typename T>
batch_acc<T> &batch_acc<T>::operator=(const batch_acc &other)
{
    size_ = other.size_;
    num_batches_ = other.num_batches_;
    base_size_ = other.base_size_;
    store_.reset(other.store_ ? new batch_data<T>(*other.store_) : nullptr);
    cursor_ = other.cursor_;
    return *this;
}

template <typename T>
void batch_acc<T>::reset()
{
    cursor_ = 0;
    if (valid())
        store_->reset();
    else
        store_.reset(new batch_data<T>(size_, num_batches_));
    store_->set_batch_size(base_size_);
}

template <typename T>
void batch_acc<T>::set_size(size_t size)
{
    size_ = size;
    if (valid()) {
        store_.reset(new batch_data<T>(size_, num_batches_));
        reset();
    }
}

template <typename T>
void batch_acc<T>::set_batch_size(uint64_t batch_size)
{
    base_size_ = batch_size;
    if (valid())
        reset();
}

template <typename T>
void batch_acc<T>::set_num_batches(size_t num_batches)
{
    // TODO: handle the case where we just discard levels more gracefully
    num_batches_ = num_batches;
    if (valid()) {
        store_.reset(new batch_data<T>(size_, num_batches_));
        reset();
    }
}

template <typename T>
void batch_acc<T>::add(const computed<T> &source, uint64_t count)
{
    internal::check_valid(*this);

    // All batches are full: add neighbours to make room in the upper half.
    if (cursor_ == num_batches_) {
        batch_list<T> batches(*store_);
        batches.coarsen();
        store_->reset();
        cursor_ = batches.store(*store_);
    }

    // Since Eigen matrix are column-major, we can just pass the pointer
    source.add_to(view<T>(store_->batch().col(cursor_).data(), size()));
    store_->count()(cursor_) += count;
    if (store_->count()(cursor_) >= store_->batch_size())
        ++cursor_;
}

template <typename T>
batch_acc<T> &batch_acc<T>::operator<<(const batch_result<T> &other)
{
    internal::check_valid(*this);
    if (size() != other.size())
        throw size_mismatch();

    batch_list<T> mine(*store_), theirs(other.store());
    const uint64_t batch_size = std::max(mine.batch_size(), theirs.batch_size());
    mine.coarsen_to(batch_size);
    theirs.coarsen_to(batch_size);
    mine.append(theirs);
    mine.fit(num_batches_);

    store_->reset();
    cursor_ = mine.store(*store_);
    return *this;
}

template <typename T>
batch_result<T> batch_acc<T>::result() const
{
    internal::check_valid(*this);
    batch_result<T> result(*store_);
    return result;
}

template <typename T>
batch_result<T> batch_acc<T>::finalize()
{
    batch_result<T> result;
    finalize_to(result);
    return result;
}

template <typename T>
void batch_acc<T>::finalize_to(batch_result<T> &result)
{
    internal::check_valid(*this);
    result.store_.reset();
    result.store_.swap(store_);
}

template class batch_acc<double>;
template class batch_acc<std::complex<double> >;


template <typename T>
batch_result<T>::batch_result(const batch_result &other)
    : store_(other.store_ ? new batch_data<T>(*other.store_) : nullptr)
{ }

template <typename T>
batch_result<T> &batch_result<T>::operator=(const batch_result &other)
{
    store_.reset(other.store_ ? new batch_data<T>(*other.store_) : nullptr);
    return *this;
}

template <typename T>
bool operator==(const batch_result<T> &r1, const batch_result<T> &r2)
{
    return r1.count() == r2.count()
        && r1.store().batch() == r2.store().batch();
}

template bool operator==(const batch_result<double> &r1,
                         const batch_result<double> &r2);
template bool operator==(const batch_result<std::complex<double>> &r1,
                         const batch_result<std::complex<double>> &r2);

template <typename T>
column<T> batch_result<T>::mean() const
{
    return store_->batch().rowwise().sum() / count();
}

template <typename T>
template <typename Str>
column<typename bind<Str,T>::var_type> batch_result<T>::var() const
{
    var_acc<T, Str> aux_acc(store_->size());
    for (size_t i = 0; i != store_->num_batches(); ++i) {
        aux_acc.add(make_adapter(store_->batch().col(i)), store_->count()(i),
                    nullptr);
    }
    return aux_acc.finalize().var();
}

template <typename T>
template <typename Str>
typename eigen<typename bind<Str,T>::cov_type>::matrix batch_result<T>::cov() const
{
    cov_acc<T, Str> aux_acc(store_->size());
    for (size_t i = 0; i != store_->num_batches(); ++i)
        aux_acc.add(make_adapter(store_->batch().col(i)), store_->count()(i));
    return aux_acc.finalize().cov();
}

template <typename T>
column<typename bind<circular_var,T>::var_type> batch_result<T>::stderror() const
{
    var_acc<T, circular_var> aux_acc(store_->size());
    for (size_t i = 0; i != store_->num_batches(); ++i) {
        aux_acc.add(make_adapter(store_->batch().col(i)), store_->count()(i),
                    nullptr);
    }
    return aux_acc.finalize().stderror();
}

template <typename T>
void batch_result<T>::reduce(const reducer &r, bool pre_commit, bool post_commit)
{
    internal::check_valid(*this);
    const reducer_setup setup = r.get_setup();

    // Each instance writes its batches to its own columns of a wider store, so
    // that the sum over instances puts them side by side: `slots` columns for
    // full batches and one for the partial batch.
    if (pre_commit) {
        reduce_num_batches_ = num_batches();
        size_t slots = reduce_num_batches_;
        if (setup.count > 1)
            slots = std::max<size_t>(2, std::min<size_t>(slots, 2 * slots / setup.count));

        batch_list<T> mine(*store_);
        mine.fit(slots);
        mine.coarsen_to(r.get_max(mine.batch_size()));

        store_.reset(new batch_data<T>(size(), setup.count * (slots + 1)));
        mine.store(*store_, setup.pos * (slots + 1));
        r.reduce(view<T>(store_->batch().data(), store_->batch().size()));
        r.reduce(view<uint64_t>(store_->count().data(), store_->num_batches()));
    }
    if (pre_commit && post_commit) {
        r.commit();
    }
    if (post_commit) {
        if (setup.have_result) {
            batch_list<T> all(*store_);
            all.fit(reduce_num_batches_);
            store_.reset(new batch_data<T>(size(), reduce_num_batches_));
            all.store(*store_);
        } else {
            store_.reset();   // free data
        }
    }
}

template column<double> batch_result<double>::var<circular_var>() const;
template column<double> batch_result<std::complex<double> >::var<circular_var>() const;
template column<complex_op<double> > batch_result<std::complex<double> >::var<elliptic_var>() const;

template eigen<double>::matrix batch_result<double>::cov< circular_var>() const;
template eigen<std::complex<double>>::matrix batch_result<std::complex<double> >::cov<circular_var>() const;
template eigen<complex_op<double> >::matrix batch_result<std::complex<double> >::cov<elliptic_var>() const;

template class batch_result<double>;
template class batch_result<std::complex<double> >;


template <typename T>
void serialize(serializer &s, const std::string &key, const batch_result<T> &self)
{
    internal::check_valid(self);
    internal::serializer_sentry group(s, key);

    // Serialize as 64-bit integers for consistency
    serialize(s, "@size", static_cast<uint64_t>(self.size()));
    serialize(s, "@num_batches", static_cast<uint64_t>(self.store().num_batches()));
    serialize(s, "@batch_size", self.store().batch_size());

    s.enter("batch");
    serialize(s, "count", self.store().count());
    serialize(s, "sum", self.store().batch());
    s.exit();

    s.enter("mean");
    serialize(s, "value", self.mean());
    serialize(s, "error", self.stderror());
    s.exit();
}

template <typename T>
void deserialize(deserializer &s, const std::string &key, batch_result<T> &self)
{
    typedef typename bind<circular_var, T>::var_type var_type;
    internal::deserializer_sentry group(s, key);

    // first deserialize the fundamentals and make sure that the target fits
    uint64_t new_size, new_nbatches, new_batch_size;
    deserialize(s, "@size", new_size);
    deserialize(s, "@num_batches", new_nbatches);
    deserialize(s, "@batch_size", new_batch_size);
    if (!self.valid() || self.size() != new_size || self.store().num_batches() != new_nbatches)
        self.store_.reset(new batch_data<T>(new_size, new_nbatches));

    self.store().set_batch_size(new_batch_size);

    // deserialize data
    s.enter("batch");
    deserialize(s, "count", self.store().count());
    deserialize(s, "sum", self.store().batch());
    s.exit();

    size_t new_size_sizet = new_size;
    s.enter("mean");
    s.read("value", ndview<T>(nullptr, &new_size_sizet, 1)); // discard
    s.read("error", ndview<var_type>(nullptr, &new_size_sizet, 1)); // discard
    s.exit();
}

template void serialize(serializer &, const std::string &key, const batch_result<double> &);
template void serialize(serializer &, const std::string &key, const batch_result<std::complex<double>> &);

template void deserialize(deserializer &, const std::string &key, batch_result<double> &);
template void deserialize(deserializer &, const std::string &key, batch_result<std::complex<double> > &);

template <typename T>
std::ostream &operator<<(std::ostream &str, const batch_result<T> &self)
{
    internal::check_valid(self);
    internal::format_sentry sentry(str);
    verbosity verb = internal::get_format(str, PRINT_TERSE);

    if (verb == PRINT_VERBOSE)
        str << "<X> = ";
    str << self.mean() << " +- " << self.stderror();

    if (verb == PRINT_VERBOSE) {
        str << "\n<Xi> = " << self.store().batch()
            << "\nNi = " << self.store().count();
    }
    return str;
}

template std::ostream &operator<<(std::ostream &, const batch_result<double> &);
template std::ostream &operator<<(std::ostream &, const batch_result<std::complex<double>> &);

}} /* namespace alps::alea */
