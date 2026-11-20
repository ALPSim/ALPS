/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include <alps/alea/batch.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/serialize.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/state.hpp>
#include <alps/alea/internal/format.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace alps { namespace alea {

namespace {

/** Batch size of a store; if not recorded, its smallest non-empty batch */
template <typename T>
uint64_t stored_batch_size(const batch_data<T> &data, size_t begin, size_t end)
{
    uint64_t batch_size = data.batch_size();
    if (batch_size == 0) {
        for (size_t i = begin; i != end; ++i) {
            if (data.count()(i) != 0 && (batch_size == 0 || data.count()(i) < batch_size))
                batch_size = data.count()(i);
        }
        if (batch_size == 0)
            batch_size = 1;
    }
    return batch_size;
}

/**
 * The batches of a series as a list: full batches in order, and the samples
 * left over at the end pooled into one partial batch. The full batches are
 * the columns of one matrix, so that coarsening works in place.
 */
template <typename T>
class batch_list
{
public:
    batch_list(size_t size, uint64_t batch_size, size_t capacity=0)
        : batch_size_(batch_size)
        , full_(size, capacity)
        , count_(capacity)
        , num_full_(0)
        , partial_(column<T>::Zero(size))
        , partial_count_(0)
    { }

    /** Read the non-empty batches in the given range of columns */
    batch_list(const batch_data<T> &data, size_t begin, size_t end, size_t capacity=0)
        : batch_list(data.size(), stored_batch_size(data, begin, end), std::max(capacity, end - begin))
    {
        for (size_t i = begin; i != end; ++i)
            add(data.batch().col(i), data.count()(i));
    }

    explicit batch_list(const batch_data<T> &data, size_t capacity=0)
        : batch_list(data, 0, data.num_batches(), capacity)
    { }

    uint64_t batch_size() const { return batch_size_; }

    size_t num_full() const { return num_full_; }

    template <typename Column>
    void add(const Column &sum, uint64_t count)
    {
        if (count >= batch_size_) {
            push(sum, count);
        } else if (count != 0) {
            partial_ += sum;
            partial_count_ += count;
            promote();
        }
    }

    void append(const batch_list &other)
    {
        reserve(num_full_ + other.num_full_ + 1);
        for (size_t i = 0; i != other.num_full_; ++i)
            add(other.full_.col(i), other.count_[i]);
        add(other.partial_, other.partial_count_);
    }

    /**
     * Append the batches of a store of the same batch size without copying
     * it first: the same as appending `batch_list(data)`.
     */
    void append(const batch_data<T> &data)
    {
        assert(stored_batch_size(data, 0, data.num_batches()) == batch_size_);
        reserve(num_full_ + data.num_batches() + 1);
        column<T> partial = column<T>::Zero(partial_.rows());
        uint64_t partial_count = 0;
        for (size_t i = 0; i != data.num_batches(); ++i) {
            const uint64_t count = data.count()(i);
            if (count >= batch_size_) {
                add(data.batch().col(i), count);
            } else if (count != 0) {
                partial += data.batch().col(i);
                partial_count += count;
                if (partial_count >= batch_size_) {
                    add(partial, partial_count);
                    partial.setZero();
                    partial_count = 0;
                }
            }
        }
        add(partial, partial_count);
    }

    /** Add neighbouring batches and double the batch size */
    void coarsen()
    {
        const size_t pairs = num_full_ / 2;
        for (size_t i = 0; i != pairs; ++i) {
            full_.col(i) = full_.col(2 * i) + full_.col(2 * i + 1);
            count_[i] = count_[2 * i] + count_[2 * i + 1];
        }
        if (num_full_ % 2 != 0) {
            partial_ += full_.col(num_full_ - 1);
            partial_count_ += count_[num_full_ - 1];
        }
        num_full_ = pairs;
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
        while (num_full_ > num_batches
               || (num_full_ == num_batches && partial_count_ != 0))
            coarsen();
    }

    /**
     * Write to the store, which must have room, and clear its other columns;
     * returns the number of full batches
     */
    size_t store(batch_data<T> &data) const
    {
        size_t used = num_full_;
        data.batch().leftCols(num_full_) = full_.leftCols(num_full_);
        for (size_t i = 0; i != num_full_; ++i)
            data.count()(i) = count_[i];
        if (partial_count_ != 0) {
            data.batch().col(used) = partial_;
            data.count()(used) = partial_count_;
            ++used;
        }
        const size_t rest = data.num_batches() - used;
        data.batch().rightCols(rest).setZero();
        data.count().tail(rest).setZero();
        data.set_batch_size(batch_size_);
        return num_full_;
    }

    /** Full batches to the first columns and the partial batch to column `last` */
    void store_aligned(batch_data<T> &data, size_t last) const
    {
        if (num_full_ > last)
            throw std::logic_error("More batches than columns");
        data.batch().leftCols(num_full_) = full_.leftCols(num_full_);
        for (size_t i = 0; i != num_full_; ++i)
            data.count()(i) = count_[i];
        data.batch().col(last) = partial_;
        data.count()(last) = partial_count_;
        data.set_batch_size(batch_size_);
    }

private:
    template <typename Column>
    void push(const Column &sum, uint64_t count)
    {
        if (num_full_ == size_t(full_.cols()))
            reserve(std::max<size_t>(2 * num_full_, 1));
        full_.col(num_full_) = sum;
        count_[num_full_] = count;
        ++num_full_;
    }

    void reserve(size_t capacity)
    {
        if (capacity <= size_t(full_.cols()))
            return;
        full_.conservativeResize(Eigen::NoChange, capacity);
        count_.resize(capacity);
    }

    void promote()
    {
        if (partial_count_ < batch_size_)
            return;
        push(partial_, partial_count_);
        partial_.setZero();
        partial_count_ = 0;
    }

    uint64_t batch_size_;
    typename eigen<T>::matrix full_;
    std::vector<uint64_t> count_;
    size_t num_full_;
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

    const batch_data<T> &theirs = other.store();
    const uint64_t their_batch_size = stored_batch_size(theirs, 0, theirs.num_batches());
    batch_list<T> mine(*store_, num_batches_ + theirs.num_batches() + 1);
    const uint64_t batch_size = std::max(mine.batch_size(), their_batch_size);
    mine.coarsen_to(batch_size);
    if (their_batch_size == batch_size) {
        mine.append(theirs);
    } else {
        batch_list<T> coarsened(theirs);
        coarsened.coarsen_to(batch_size);
        mine.append(coarsened);
    }
    mine.fit(num_batches_);
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

template <typename T>
void batch_acc<T>::save_state(serializer &s, const std::string &key) const
{
    internal::check_valid(*this);
    internal::serializer_sentry group(s, key);
    internal::write_state_header(s, internal::STATE_BATCH);
    serialize(s, "@size", static_cast<uint64_t>(size_));
    serialize(s, "@num_batches", static_cast<uint64_t>(num_batches_));
    serialize(s, "@base_size", static_cast<uint64_t>(base_size_));
    serialize(s, "@batch_size", store_->batch_size());
    serialize(s, "cursor", static_cast<uint64_t>(cursor_));
    serialize(s, "count", store_->count());
    serialize(s, "sum", store_->batch());
}

template <typename T>
void batch_acc<T>::load_state(deserializer &s, const std::string &key)
{
    internal::deserializer_sentry group(s, key);
    internal::read_state_header(s, internal::STATE_BATCH);
    size_ = internal::read_uint(s, "@size");
    num_batches_ = internal::read_uint(s, "@num_batches");
    base_size_ = internal::read_uint(s, "@base_size");
    store_.reset(new batch_data<T>(size_, num_batches_));
    store_->set_batch_size(internal::read_uint(s, "@batch_size"));
    cursor_ = internal::read_uint(s, "cursor");
    deserialize(s, "count", store_->count());
    deserialize(s, "sum", store_->batch());
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

    // The runs are independent, so the i-th batches of all instances add up to
    // one batch of the combined run. All instances coarsen to a common batch
    // size and the sums are taken column by column: the buffer has the same
    // `num_batches() + 1` columns at any number of instances (the extra column
    // pools the partial batches).
    if (pre_commit) {
        reduce_num_batches_ = r.get_max(num_batches());
        batch_list<T> mine(*store_);
        mine.fit(reduce_num_batches_);
        mine.coarsen_to(r.get_max(mine.batch_size()));

        store_.reset(new batch_data<T>(size(), reduce_num_batches_ + 1));
        mine.store_aligned(*store_, reduce_num_batches_);
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
