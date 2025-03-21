/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Ilia Zintchenko <iliazin@gmail.com>                       *
 *                       Jan Gukelberger                                           *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_MULTI_ARRAY_BASE_HPP
#define ALPS_MULTI_ARRAY_BASE_HPP

#include <boost/multi_array.hpp>

namespace alps{

    template <class T, std::size_t D, class Allocator = std::allocator<T> >
    class multi_array : public boost::multi_array<T,D,Allocator>
    {
        typedef boost::multi_array<T,D,Allocator> base_type;
        typedef multi_array<T,D,Allocator> array_type;

    public:
        typedef typename base_type::size_type size_type;

        multi_array(size_type N, size_type M, size_type K, size_type J, size_type I) : base_type(boost::extents[N][M][K][J][I]) {}
        multi_array(size_type N, size_type M, size_type K, size_type J) : base_type(boost::extents[N][M][K][J]) {}
        multi_array(size_type N, size_type M, size_type K) : base_type(boost::extents[N][M][K]) {}
        multi_array(size_type N, size_type M) : base_type(boost::extents[N][M]) {}
        multi_array(size_type N) : base_type(boost::extents[N]) {}


	multi_array(multi_array const& a) : base_type(static_cast<base_type const&>(a)) {}
	multi_array(base_type const& a) : base_type(a) {}
    multi_array(const boost::detail::multi_array::extent_gen<D>& ext) : base_type(ext) {}


        multi_array() : base_type() {}

        template<typename Array>
        multi_array(const Array& x) : base_type(x) {};


        array_type& operator=(const array_type& a)
        {
            if(this != &a)
            {
                const size_type* shp = a.shape();
                std::vector<size_type> ext(shp,shp+a.num_dimensions());
                (*this).resize(ext);
                base_type::operator=(a);
            }

            return *this;
        }

        template <typename Array>
        multi_array& operator=(const Array& a) 
        {
            const typename Array::size_type* shp = a.shape();
            std::vector<size_type> ext(shp,shp+a.num_dimensions());
            (*this).resize(ext);
            base_type::operator=(a);
            return *this;
        }

        array_type& operator+=(const array_type& a)
        {
            // TODO: how do we handle that?
            if (std::accumulate(this->shape(), this->shape() + D, size_type(0)) == 0) {
                boost::array<T, D> extent;
                std::copy(a.shape(), a.shape() + D, extent.begin());
                this->resize(extent);
            }
            assert(std::equal(this->shape(),this->shape()+D,a.shape()));
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),a.data(),(*this).data(),std::plus<T>());
            return *this;
        }

        array_type& operator+=(const T s)
        {
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),(*this).data(),std::bind2nd(std::plus<T>(),s));
            return *this;
        }

        array_type& operator-=(const array_type& a)
        {
            assert(std::equal(this->shape(),this->shape()+D,a.shape()));
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),a.data(),(*this).data(),std::minus<T>());
            return *this;
        }

        array_type& operator-=(const T s)
        {
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),(*this).data(),std::bind2nd(std::minus<T>(),s));
            return *this;
        }

        array_type& operator*=(const array_type& a)
        {
            assert(std::equal(this->shape(),this->shape()+D,a.shape()));
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),a.data(),(*this).data(),std::multiplies<T>());
            return *this;
        }

        array_type& operator*=(const T s)
        {
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),(*this).data(),std::bind2nd(std::multiplies<T>(),s));
            return *this;
        }

        array_type& operator/=(const array_type& a)
        {
            assert(std::equal(this->shape(),this->shape()+D,a.shape()));
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),a.data(),(*this).data(),std::divides<T>());
            return *this;
        }

        array_type& operator/=(const T s)
        {
            std::transform((*this).data(),(*this).data()+(*this).num_elements(),(*this).data(),std::bind2nd(std::divides<T>(),s));
            return *this;
        }

    };//class multi_array

}//namespace alps

#endif // ALPS_MULTI_ARRAY_BASE_HPP
