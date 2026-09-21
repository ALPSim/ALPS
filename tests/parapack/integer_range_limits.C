// SPDX-License-Identifier: MIT

#include <alps/parapack/integer_range.h>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

template <class T>
void check_limits() {
  using range = alps::integer_range<T>;
  const T low = (std::numeric_limits<T>::min)();
  const T high = (std::numeric_limits<T>::max)();
  constexpr int bits = std::numeric_limits<T>::digits + std::numeric_limits<T>::is_signed;
  const range empty;
  require(empty.empty() && !empty.valid() && empty.size() == 0, "Default range is not empty");
  require(range(low).size() == 1 && range(high).size() == 1, "Limit singleton has wrong size");
  range full(low, high);
  require(full.valid() && !full.empty(), "Full range is incorrectly empty");
  full.include(T(0));
  require(full.min() == low && full.max() == high, "Including zero collapsed a full range");
  require(full.overlap(range(T(0), T(3))).size() == 4, "Full range overlap failed");
  if constexpr (bits < std::numeric_limits<std::uintmax_t>::digits) {
    require(full.size() == (std::uintmax_t(1) << bits), "Full range has wrong size");
  } else {
    bool overflow = false;
    try {
      (void)full.size();
    } catch (const std::overflow_error&) {
      overflow = true;
    }
    require(overflow, "Unrepresentable size did not report overflow");
    require(range(low, high - 1).size() == (std::numeric_limits<std::uintmax_t>::max)(),
            "Largest representable size is incorrect");
  }
}

int main() {
  try {
    check_limits<std::int32_t>();
    check_limits<std::uint32_t>();
    check_limits<std::intmax_t>();
    check_limits<std::uintmax_t>();
    require(alps::integer_range<int>(-7, 9).size() == 17, "Range crossing zero has wrong size");
    require(alps::integer_range<int>(-9, -7).size() == 3, "Negative range has wrong size");
    require(alps::integer_range<int>(3, 1).empty(), "Reversed range is not empty");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
