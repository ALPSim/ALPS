// SPDX-License-Identifier: MIT

// This consumer must compile with the utilities/container include roots only.
// It exercises array storage and formatting without the numerical umbrella.
#include <alps/multi_array/multi_array.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/utility/resize.hpp>

#include <iostream>
#include <sstream>

int main()
{
    alps::multi_array<int, 2> values(2, 3);
    for (std::size_t i = 0; i < values.num_elements(); ++i)
        values.data()[i] = static_cast<int>(i + 1);

    alps::multi_array<int, 2> copy;
    copy = values;
    if (copy.shape()[0] != 2 || copy.shape()[1] != 3 || copy[1][2] != 6) {
        std::cerr << "array assignment did not preserve shape and values\n";
        return 1;
    }
    copy[0][0] = 9;
    if (values[0][0] != 1) {
        std::cerr << "array assignment did not create independent storage\n";
        return 1;
    }

    alps::multi_array<double, 2> resized(1, 1);
    alps::resize_same_as(resized, values);
    if (resized.shape()[0] != 2 || resized.shape()[1] != 3) {
        std::cerr << "array resize did not preserve the source shape\n";
        return 1;
    }

    std::ostringstream output;
    output << alps::short_print(values);
    if (output.str() != "[1,..6..,6]") {
        std::cerr << "array short_print changed: " << output.str() << '\n';
        return 1;
    }
}
