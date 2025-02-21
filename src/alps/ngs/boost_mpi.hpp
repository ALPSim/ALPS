/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_BOOST_MPI_HPP
#define ALPS_NGS_BOOST_MPI_HPP

#ifdef ALPS_HAVE_MPI

    #include <boost/mpi.hpp>
    #include <boost/array.hpp>
    #include <alps/multi_array.hpp>

    #include <vector>

    namespace boost {
        namespace mpi {

            // std::vectyor
            template<typename T, typename A, typename Op> void reduce(const communicator & comm, std::vector<T, A> const & in_values, Op op, int root) {
                reduce(comm, &in_values.front(), in_values.size(), op, root);
            }

            template<typename T, typename A, typename Op> void reduce(const communicator & comm, std::vector<T, A> const & in_values, std::vector<T, A> & out_values, Op op, int root) {
                out_values.resize(in_values.size());
                reduce(comm, &in_values.front(), in_values.size(), &out_values.front(), op, root);
            }

            // boost::array
            template<typename T, std::size_t N, typename Op> void reduce(const communicator & comm, boost::array<T, N> const & in_values, Op op, int root) {
                reduce(comm, &in_values.front(), in_values.size(), op, root);
            }

            template<typename T, std::size_t N, typename Op> void reduce(const communicator & comm, boost::array<T, N> const & in_values, boost::array<T, N> & out_values, Op op, int root) {
                reduce(comm, &in_values.front(), in_values.size(), &out_values.front(), op, root);
            }

            // boost::mulit_array
            template<typename T, std::size_t N, typename A, typename Op> void reduce(const communicator & comm, boost::multi_array<T, N, A> const & in_values, Op op, int root) {
                reduce(comm, in_values.data(), std::accumulate(in_values.shape(), in_values.shape() + boost::multi_array<T, N, A>::dimensionality, 0), op, root);
            }

            template<typename T, std::size_t N, typename A, typename Op> void reduce(const communicator & comm, boost::multi_array<T, N, A> const & in_values, boost::multi_array<T, N, A> & out_values, Op op, int root) {
                boost::array<typename boost::multi_array<T, N, A>::size_type, boost::multi_array<T, N, A>::dimensionality> shape;
                std::copy(in_values.shape(), in_values.shape() + boost::multi_array<T, N, A>::dimensionality, shape.begin());
                out_values.resize(shape);
                reduce(comm, in_values.data(), std::accumulate(shape.begin(), shape.end(), 0), out_values.data(), op, root);
            }

            // alps::mulit_array
            template<typename T, std::size_t N, typename A, typename Op> void reduce(const communicator & comm, alps::multi_array<T, N, A> const & in_values, Op op, int root) {
                reduce(comm, static_cast<boost::multi_array<T, N, A> const &>(in_values), op, root);
            }

            template<typename T, std::size_t N, typename A, typename Op> void reduce(const communicator & comm, alps::multi_array<T, N, A> const & in_values, alps::multi_array<T, N, A> & out_values, Op op, int root) {
                reduce(comm, static_cast<boost::multi_array<T, N, A> const &>(in_values), static_cast<boost::multi_array<T, N, A> &>(out_values), op, root);
            }

        }
    }

#endif

#endif
