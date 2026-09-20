// SPDX-License-Identifier: MIT
#include "tinyvector/tinyvector.hpp"
#include <cstdio>
#include <stdexcept>

int main() {
    const tinyvector<double, 3> a(std::vector<double>{1, 2, 3});
    const tinyvector<double, 3> b(std::vector<double>{4, -5, 6});
    if (dot(a, b) != 12 || dot(a, a) != 14 || dot(a - a, b) != 0)
        throw std::runtime_error("incorrect spin dot product");
    auto scaled = (a + b) / 2.0;
    if (scaled[0] != 2.5 || scaled[1] != -1.5 || scaled[2] != 4.5)
        throw std::runtime_error("incorrect spin arithmetic");
    if (std::vector<double>(a.begin(), a.end()) != tinyvector<double, 3>::vector(a))
        throw std::runtime_error("incorrect const spin iteration");
    {
        alps::hdf5::archive output("tinyvector-test.h5", "w");
        output["spin"] << scaled;
    }
    tinyvector<double, 3> restored;
    {
        alps::hdf5::archive input("tinyvector-test.h5", "r");
        input["spin"] >> restored;
    }
    std::remove("tinyvector-test.h5");
    if (dot(restored - scaled, restored - scaled) != 0)
        throw std::runtime_error("incorrect spin checkpoint round trip");
}
