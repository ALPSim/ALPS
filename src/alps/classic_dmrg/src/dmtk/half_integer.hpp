// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The half-integer value type of dmtk's quantum numbers: the subset of
// <alps/model/half_integer.h> (this tree) that dmtk uses, without the Boost
// type-trait machinery, so that the classic DMRG component does not include
// the legacy model headers.
#pragma once

#include <cmath>
#include <ostream>
#include <type_traits>

namespace dmtk {

template<class I>
class half_integer {
public:
    using integer_type = I;
    half_integer() = default;
    template<class J, class = std::enable_if_t<std::is_integral<J>::value>>
    half_integer(J x) : val_(static_cast<I>(2 * x)) {}
    half_integer(double x) : val_(static_cast<I>(std::floor(2 * x + 0.5))) {}

    double to_double() const { return 0.5 * val_; }
    I get_twice() const { return val_; }
    void set_half(I twice) { val_ = twice; }

    bool operator==(half_integer const& o) const { return val_ == o.val_; }
    bool operator!=(half_integer const& o) const { return val_ != o.val_; }
    bool operator<(half_integer const& o) const { return val_ < o.val_; }
    bool operator>(half_integer const& o) const { return val_ > o.val_; }
    bool operator<=(half_integer const& o) const { return val_ <= o.val_; }
    bool operator>=(half_integer const& o) const { return val_ >= o.val_; }

    half_integer operator-() const { half_integer r; r.val_ = static_cast<I>(-val_); return r; }
    half_integer& operator+=(half_integer const& o) { val_ = static_cast<I>(val_ + o.val_); return *this; }
    half_integer& operator-=(half_integer const& o) { val_ = static_cast<I>(val_ - o.val_); return *this; }
    half_integer operator+(half_integer const& o) const { half_integer r(*this); r += o; return r; }
    half_integer operator-(half_integer const& o) const { half_integer r(*this); r -= o; return r; }

private:
    I val_ = 0;
};

template<class I>
std::ostream& operator<<(std::ostream& out, half_integer<I> const& x) {
    if (x.get_twice() % 2 == 0) return out << x.get_twice() / 2;
    return out << x.get_twice() << "/2";
}

} // namespace dmtk
