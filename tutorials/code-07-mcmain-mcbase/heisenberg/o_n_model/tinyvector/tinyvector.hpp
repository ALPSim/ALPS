/* tinyvector.hpp
 *
 * contains the default implementation
 */

#ifndef TINYVECTOR_HPP
#define TINYVECTOR_HPP

#include <boost/array.hpp>
#include <alps/hdf5.hpp>
#include <alps/hdf5/array.hpp>

#include <vector>
#include <cmath>
#include <iostream>

template <class T = double, int N = 3>
class tinyvector {
    public:
        typedef typename boost::array<T, N> data_type;
        typedef typename data_type::value_type value_type;
        typedef typename data_type::iterator iterator;
        typedef typename data_type::const_iterator const_iterator;

        tinyvector() : _data() {}
        tinyvector(T value) : _data() { initialize(value); }
        tinyvector(const std::vector<T> &vec) {
            for(int i = 0; i < N; ++i)
                _data[i] = vec[i];
        }
        const value_type * data() const { return _data.data(); }
        value_type * data() { return _data.c_array(); }

        const value_type & front() const { return _data.front(); }
        value_type & front() { return _data.front(); }

        const_iterator begin() const { return _data.begin(); }
        iterator begin() { return _data.begin(); }
        const_iterator end() const { return _data.end(); }
        iterator end() { return _data.end(); }

        inline const T operator[](int i) const {
            return _data[i]; }
        inline T & operator[](int i) {
            return _data[i]; }

        void initialize(T init) {
            for(int i = 0; i < N; ++i)
                _data[i] = init;
        }

        static inline const std::vector<T> vector(const tinyvector<T, N> & tv){
            std::vector<T> result;
            for(int i = 0; i < N; ++i)
                result.push_back(tv[i]);
            return result;
        }

        void save(alps::hdf5::archive & ar) const {
            ar << alps::make_pvp("data", _data);
        }
        void load(alps::hdf5::archive & ar) {
            ar >> alps::make_pvp("data", _data);
        }
    private:
        data_type _data;
};

template <class T, int N>
inline const tinyvector<T, N> & operator+=(tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] += right[i];
    return left;
};

template <class T, int N>
inline const tinyvector<T, N> & operator-=(tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] -= right[i];
    return left;
}

template <class T, int N>
inline const tinyvector<T, N> & operator*=(tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right[i];
    return left;
}

template <class T, int N>
inline const tinyvector<T, N> & operator*=(tinyvector<T, N> &left, T right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right;
    return left;
}

template <class T, int N>
inline const tinyvector<T, N> & operator/=(tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right[i];
    return left;
}

template <class T, int N>
inline const tinyvector<T, N> & operator/=(tinyvector<T, N> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right;
    return left;
}

template <class T, int N>
inline const tinyvector<T, N> operator+(const tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    tinyvector<T, N> result(left);
    return result += right;
}

template <class T, int N>
inline const tinyvector<T, N> operator-(const tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    tinyvector<T, N> result(left);
    return result -= right;
}

template <class T, int N>
inline const tinyvector<T, N> operator-(const tinyvector<T, N> &spin) {
    tinyvector<T, N> result;
    for(int i = 0; i < N; ++i)
        result[i] = - spin[i];
    return result;
}

template <class T, int N>
inline const tinyvector<T, N> operator*(const tinyvector<T, N> &left, double right) {
    tinyvector<T, N> result(left);
    return result *= right;
}

template <class T, int N>
inline const tinyvector<T, N> operator*(double left, const tinyvector<T, N> &right) {
    return right * left;
}

template <class T, int N>
inline const tinyvector<T, N> operator/(const tinyvector<T, N> &left, double right) {
    tinyvector<T, N> result(left);
    return result /= right;
}

template <class T, int N>
inline double sum(const tinyvector<T, N> &spin) {
    double result = 0.;
    for(int i = 0; i < N; ++i)
        result += spin[i];
    return result;
}

template <class T, int N>
inline double dot(const tinyvector<T, N> &left, const tinyvector<T, N> &right) {
    double result = 0.;
    for (int i = 0; i < N; ++i)
        result += left[i] * right[i];
    return result;
}

template <class T, int N>
inline double abs(const tinyvector<T, N> &spin) {
   return std::sqrt(dot(spin, spin));
}

template <class T, int N>
std::ostream& operator<<(std::ostream& os, const tinyvector<T, N>& spin)
{
    //~ os << "[ ";
    for(int i = 0; i < N; ++i)
        os << spin[i] << " ";
    //~ os << " ]";
    return os;
}

#endif
