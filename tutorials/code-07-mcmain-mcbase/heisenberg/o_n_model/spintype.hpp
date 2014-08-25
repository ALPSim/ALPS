/* helper.hpp for ndim_spin example
 */

#ifndef HEISENBERG_HELPER_HPP
#define HEISENBERG_HELPER_HPP

//~ #include "ndim_spin.hpp"

#include <boost/array.hpp>
#include <alps/hdf5.hpp>

#include <vector>
#include <cmath>
#include <iostream>

template <int N>
class spintype {
    public:
        inline const double & operator[](int i) const { return _data[i]; };
        inline double & operator[](int i) { return _data[i]; };
        void initialize(double init) {
            for(int i = 0; i < N; ++i)
                _data[i] = init;
        }
        void save(alps::hdf5::archive & ar) const {
            ar << alps::make_pvp("data", _data);
        }
        void load(alps::hdf5::archive & ar) {
            ar >> alps::make_pvp("data", _data);
        }
    private:
        boost::array<double, N> _data;
};

template <int N>
const spintype<N> & operator+=(spintype<N> &left, const spintype<N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] += right[i];
    return left;
}

template <int N>
const spintype<N> & operator-=(spintype<N> &left, const spintype<N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] -= right[i];
    return left;
}

template <int N>
const spintype<N> & operator*=(spintype<N> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right;
    return left;
}

template <int N>
const spintype<N> & operator/=(spintype<N> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right;
    return left;
}

template <int N>
const spintype<N> operator+(const spintype<N> &left, const spintype<N> &right) {
    spintype<N> result(left);
    return result += right;
}

template <int N>
const spintype<N> operator-(const spintype<N> &left, const spintype<N> &right) {
    spintype<N> result(left);
    return result -= right;
}

template <int N>
const spintype<N> operator-(const spintype<N> &spin) {
    spintype<N> result;
    for(int i = 0; i < N; ++i)
        result[i] = - spin[i];
    return result;
}

template <int N>
const spintype<N> operator*(const spintype<N> &left, double right) {
    spintype<N> result(left);
    return result *= right;
}

template <int N>
const spintype<N> operator*(double left, const spintype<N> &right) {
    return right * left;
}

template <int N>
const spintype<N> operator/(const spintype<N> &left, double right) {
    spintype<N> result(left);
    return result /= right;
}

template <int N>
double dot(const spintype<N> &left, const spintype<N> &right) {
    double result = 0.;
    for(int i = 0; i < N; ++i)
        result += left[i] * right[i];
    return result;
}

template <int N>
double abs(const spintype<N> &spin) {
   return std::sqrt(dot(spin, spin));
}

template <int N>
std::ostream& operator<<(std::ostream& os, const spintype<N>& spin)
{
    os << "[ ";
    for(int i = 0; i < N; ++i)
        os << spin[i] << " ";
    os << " ]";
    return os;
}

template <int N>
const std::vector<double> vector_from_spintype(const spintype<N> &spin) {
    std::vector<double> result;
    for(int i = 0; i < N; ++i)
        result.push_back(spin[i]);
    return result;
}

template <int N>
const spintype<N> spin_from_vector(const std::vector<double> &vec) {
    spintype<N> result;
    for(int i = 0; i < N; ++i)
        result[i] = vec[i];
    return result;
}

#endif
