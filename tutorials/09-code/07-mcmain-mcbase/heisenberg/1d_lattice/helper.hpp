/* helper.hpp for heisenberg example
 */

#ifndef HEISENBERG_HELPER_HPP
#define HEISENBERG_HELPER_HPP

#include "heisenberg.hpp"

#include <boost/array.hpp>

#include <vector>
#include <cmath>
#include <iostream>

typedef boost::array<double, 3> spintype;

const spintype & operator+=(spintype &left, const spintype &right) {
    for(int i = 0; i < 3; ++i)
        left[i] += right[i];
    return left;
}

const spintype & operator-=(spintype &left, const spintype &right) {
    for(int i = 0; i < 3; ++i)
        left[i] -= right[i];
    return left;
}

const spintype & operator*=(spintype &left, double right) {
    for(int i = 0; i < 3; ++i)
        left[i] *= right;
    return left;
}

const spintype & operator/=(spintype &left, double right) {
    for(int i = 0; i < 3; ++i)
        left[i] /= right;
    return left;
}

const spintype operator+(const spintype &left, const spintype &right) {
    spintype result(left);
    return result += right;
}

const spintype operator-(const spintype &left, const spintype &right) {
    spintype result(left);
    return result -= right;
}

const spintype operator-(const spintype &spin) {
    spintype result = { {0, 0, 0} };
    return result -= spin;
}

const spintype operator*(const spintype &left, double right) {
    spintype result(left);
    return result *= right;
}

const spintype operator*(double left, const spintype &right) {
    return right * left;
}

const spintype operator/(const spintype &left, double right) {
    spintype result(left);
    return result /= right;
}

double dot(const spintype &left, const spintype &right) {
    double result = 0.;
    for( int i = 0; i < 3; i += 1)
    {
        result += left[i] * right[i];
    }
    return result;
}

double abs(const spintype &spin) {
   return std::sqrt(dot(spin, spin));
}

std::ostream& operator<<(std::ostream& os, const spintype& spin)
{
    os << "[ " << spin[0] << " " << spin[1] << " " << spin[2] << " ]";
    return os;
}

const std::vector<double> vector_from_spintype(const spintype &spin) {
    std::vector<double> result;
    for( int i = 0; i < 3; i += 1)
    {
        result.push_back(spin[i]);
    }
    return result;
}

const spintype spin_from_vector(const std::vector<double> &vec) {
    spintype result = { { vec[0], vec[1], vec[2] } };
    return result;
}

#endif
