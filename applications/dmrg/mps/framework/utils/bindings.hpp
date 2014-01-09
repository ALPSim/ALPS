#ifndef MATRIX_BINDINGS_H
#define MATRIX_BINDINGS_H

#include <alps/numeric/real.hpp> 

namespace maquis { namespace bindings {

    template <typename O, typename I> struct binding { 
        static O convert(const I& m){ return static_cast<O>(m); }
    };

    template<typename O, typename I> O matrix_cast(I const& input){
       return binding<O,I>::convert(input);
    }

    template <typename T>
    struct binding< std::vector<T>, alps::numeric::diagonal_matrix<T> > {
        static std::vector<T> convert(const alps::numeric::diagonal_matrix<T>& m){
            return m.get_values();
        }
    };

    template <typename T, typename S, template<class M, class SS> class C>
    struct binding< std::vector< std::vector<T> >, C<alps::numeric::diagonal_matrix<T>, S> > {
        static std::vector< std::vector<T> > convert(const C<alps::numeric::diagonal_matrix<T>, S>& m){
            std::vector< std::vector<T> > set;
            for(size_t k = 0; k < m.n_blocks(); ++k){
                set.push_back(m[k].get_values());
            }
            return set;
        }
    };

} }

namespace maquis {
    template <class T> 
    inline typename alps::numeric::real_type<T>::type real(T f){
        return alps::numeric::real(f);
    }

    template<typename _InputIterator, typename _Tp>
    inline _Tp
    accumulate(_InputIterator __first, _InputIterator __last, _Tp __init){
        return std::accumulate(__first, __last, __init);
    }

    template<typename T>
    T sqrt(T arg){
        return std::sqrt(arg);
    }
}

#endif
