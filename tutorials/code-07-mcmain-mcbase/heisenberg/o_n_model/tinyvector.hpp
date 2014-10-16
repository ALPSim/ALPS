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
class tinyvector {
    public:
        typedef typename boost::array<double, N> data_type;
        typedef typename data_type::value_type value_type;
        typedef typename data_type::iterator iterator;
        typedef typename data_type::const_iterator const_iterator;

        tinyvector() : _data() {}
        tinyvector(int dim) : _data() { assert(dim = N); }
        tinyvector(const std::vector<double> &vec) {
            for(int i = 0; i < N; ++i)
                _data[i] = vec[i];
        }

        const value_type * data() const { return _data.data(); }
        value_type * data() { return _data.c_array(); }

        const value_type & front() const { return _data.front(); }
        value_type & front() { return _data.front(); }

        const iterator begin() const { return _data.begin(); }
        iterator begin() { return _data.begin(); }
        const iterator end() const { return _data.end(); }
        iterator end() { return _data.end(); }

        inline const double operator[](int i) const { return _data[i]; };
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
        boost::array<double, N> _data __attribute__((aligned(16 * sizeof(double))));
};

template <int N>
inline const tinyvector<N> & operator+=(tinyvector<N> &left, const tinyvector<N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] += right[i];
    return left;
}

template <int N>
inline const tinyvector<N> & operator-=(tinyvector<N> &left, const tinyvector<N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] -= right[i];
    return left;
}

template <int N>
inline const tinyvector<N> & operator*=(tinyvector<N> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right;
    return left;
}

template <int N>
inline const tinyvector<N> & operator/=(tinyvector<N> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right;
    return left;
}

template <int N>
inline const tinyvector<N> operator+(const tinyvector<N> &left, const tinyvector<N> &right) {
    tinyvector<N> result(left);
    return result += right;
}

template <int N>
inline const tinyvector<N> operator-(const tinyvector<N> &left, const tinyvector<N> &right) {
    tinyvector<N> result(left);
    return result -= right;
}

template <int N>
inline const tinyvector<N> operator-(const tinyvector<N> &spin) {
    tinyvector<N> result;
    for(int i = 0; i < N; ++i)
        result[i] = - spin[i];
    return result;
}

template <int N>
inline const tinyvector<N> operator*(const tinyvector<N> &left, double right) {
    tinyvector<N> result(left);
    return result *= right;
}

template <int N>
inline const tinyvector<N> operator*(double left, const tinyvector<N> &right) {
    return right * left;
}

template <int N>
inline const tinyvector<N> operator/(const tinyvector<N> &left, double right) {
    tinyvector<N> result(left);
    return result /= right;
}

template <int N>
double dot(const tinyvector<N> &left, const tinyvector<N> &right) {
    double result = 0.;
    for(int i = 0; i < N; ++i)
        result += left[i] * right[i];
    return result;
}

template <int N>
double abs(const tinyvector<N> &spin) {
   return std::sqrt(dot(spin, spin));
}

template <int N>
std::ostream& operator<<(std::ostream& os, const tinyvector<N>& spin)
{
    os << "[ ";
    for(int i = 0; i < N; ++i)
        os << spin[i] << " ";
    os << " ]";
    return os;
}

template <int N>
inline const std::vector<double> vector_from_tinyvector(const tinyvector<N> &spin) {
    std::vector<double> result;
    for(int i = 0; i < N; ++i)
        result.push_back(spin[i]);
    return result;
}

template <int N>
inline const tinyvector<N> spin_from_vector(const std::vector<double> &vec) {
    tinyvector<N> result;
    for(int i = 0; i < N; ++i)
        result[i] = vec[i];
    return result;
}

#endif
