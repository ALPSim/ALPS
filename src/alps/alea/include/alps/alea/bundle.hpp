/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/util.hpp>

namespace alps { namespace alea {

/**
 * Bundle of measurements.
 *
 * Some accumulator, particularly the variance, are more efficient working in
 * "bundled" mode, where n measurement are bundled or batched together.
 */
template <typename T>
class bundle
{
public:
    // The sum is allocated when it is first written: an accumulator that
    // adds each sample directly never needs it, and for a large vector it is
    // as large as the vector.
    bundle(size_t size, uint64_t target) : size_(size), target_(target), count_(0) { }

    /** Clear all accumulated data */
    void reset() { if (allocated()) sum_.fill(0); count_ = 0; }

    bool is_full() { return count_ >= target_; }

    /** Number of components of the random vector (e.g., size of mean) */
    size_t size() const { return size_; }

    /** Whether the sum has been allocated; it is zero while it has not */
    bool allocated() const { return size_t(sum_.rows()) == size_; }

    uint64_t &target() { return target_; }

    const uint64_t &target() const { return target_; }

    uint64_t &count() { return count_; }

    uint64_t count() const { return count_; }

    column<T> &sum()
    {
        if (!allocated())
            sum_ = column<T>::Zero(size_);
        return sum_;
    }

    /** Empty while not allocated */
    const column<T> &sum() const { return sum_; }

private:
    size_t size_;
    column<T> sum_;
    uint64_t target_, count_;
};

}}
