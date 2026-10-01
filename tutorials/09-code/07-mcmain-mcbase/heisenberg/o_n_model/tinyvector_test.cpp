// Verify the vector operations used by the Heisenberg energy update.
#include "tinyvector/tinyvector.hpp"
#include <cmath>
#include <vector>

template <class Optimization>
bool check() {
    using vector = tinyvector<double, 3, Optimization>;
    vector left(std::vector<double>{1., 2., 3.});
    vector right(std::vector<double>{4., 5., 6.});
    vector zero(0.);
    return dot(left, right) == 32. && dot(right, left) == 32.
        && dot(left, left) == 14. && dot(left, zero) == 0.
        && dot(right - left, left) == 18.
        && std::abs(abs(left) - std::sqrt(14.)) < 1e-12;
}

int main() {
    return check<NO_OPT>() && check<INTRIN_OPT>() ? 0 : 1;
}
