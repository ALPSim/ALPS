/* tinyvector_default.hpp
 *
 * contains the default implementation
 */

#ifndef TINYVECTOR_DEFAULT_HPP
#define TINYVECTOR_DEFAULT_HPP

#include <boost/array.hpp>
#include <alps/hdf5.hpp>

#include <vector>
#include <cmath>
#include <iostream>

struct NO_OPT {};

template <int N>
struct _next_lower_power_of_2 {
    static const int value = 1 << ( N >> 1);
};

template <class T = double, int N = 3, class Opt = NO_OPT>
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
        template<class Opt2>
        tinyvector(const tinyvector<T, N, Opt2> & vec) {
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

        inline const T operator[](int i) const { 
            return _data[i]; }
        inline T & operator[](int i) { 
            return _data[i]; }

        void initialize(T init) {
            for(int i = 0; i < N; ++i)
                _data[i] = init;
        }

        static inline const std::vector<T> vector(const tinyvector<T, N, Opt> & tv){
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
        data_type _data __attribute__((aligned( 16 * sizeof(T) )));
};

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator+=(tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    for(int i = 0; i < N; ++i)
        left[i] += right[i];
    return left;
};

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator-=(tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    for(int i = 0; i < N; ++i)
        left[i] -= right[i];
    return left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator*=(tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right[i];
    return left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator*=(tinyvector<T, N, Opt> &left, T right) {
    for(int i = 0; i < N; ++i)
        left[i] *= right;
    return left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator/=(tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right[i];
    return left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> & operator/=(tinyvector<T, N, Opt> &left, double right) {
    for(int i = 0; i < N; ++i)
        left[i] /= right;
    return left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator+(const tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    tinyvector<T, N, Opt> result(left);
    return result += right;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator-(const tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    tinyvector<T, N, Opt> result(left);
    return result -= right;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator-(const tinyvector<T, N, Opt> &spin) {
    tinyvector<T, N, Opt> result;
    for(int i = 0; i < N; ++i)
        result[i] = - spin[i];
    return result;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator*(const tinyvector<T, N, Opt> &left, double right) {
    tinyvector<T, N, Opt> result(left);
    return result *= right;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator*(double left, const tinyvector<T, N, Opt> &right) {
    return right * left;
}

template <class T, int N, class Opt>
inline const tinyvector<T, N, Opt> operator/(const tinyvector<T, N, Opt> &left, double right) {
    tinyvector<T, N, Opt> result(left);
    return result /= right;
}

template <class T, int N, class Opt>
inline double sum(const tinyvector<T, N, Opt> &spin) {
    double result = 0.;
    for(int i = 0; i < N; ++i)
        result += spin[i];
    return result;
}

template <class T, int N, class Opt>
inline double dot(const tinyvector<T, N, Opt> &left, const tinyvector<T, N, Opt> &right) {
    double result = 0.;
    result = sum(left * left);
    return result;
}

template <class T, int N, class Opt>
inline double abs(const tinyvector<T, N, Opt> &spin) {
   return std::sqrt(dot(spin, spin));
}

template <class T, int N, class Opt>
std::ostream& operator<<(std::ostream& os, const tinyvector<T, N, Opt>& spin)
{
    //~ os << "[ ";
    for(int i = 0; i < N; ++i)
        os << spin[i] << " ";
    //~ os << " ]";
    return os;
}

#endif
