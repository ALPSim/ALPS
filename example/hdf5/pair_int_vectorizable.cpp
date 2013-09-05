/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/integer.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <iostream>

namespace alps {
    namespace hdf5 {

        template<> struct scalar_type<std::pair<int, int> > {
            typedef int type;
        };

        template<> struct is_continuous<std::pair<int, int> > : public boost::true_type {};

        template<> struct has_complex_elements<std::pair<int, int> > : public boost::false_type {};

        namespace detail {

            template<> struct get_extent<std::pair<int, int> > {
                static std::vector<std::size_t> apply(std::pair<int, int> const & value) {
                    return std::vector<std::size_t>(1, 2);
                }
            };

            template<> struct set_extent<std::pair<int, int> > {
                static void apply(std::pair<int, int> &, std::vector<std::size_t> const &) {}
            };

        
            template<> struct is_vectorizable<std::pair<int, int> > {
                static bool apply(std::pair<int, int> const & value) {
                    return true;
                }
            };

            template<> struct get_pointer<std::pair<int, int> > {
                static scalar_type<std::pair<int, int> >::type * apply(std::pair<int, int> & value) {
                    return &value.first;
                }
            };
        
            template<> struct get_pointer<std::pair<int, int> const> {
                static scalar_type<std::pair<int, int> >::type const * apply(std::pair<int, int> const & value) {
                    return &value.first;
                }
            };
        }

        void save(
              archive & ar
            , std::string const & path
            , std::pair<int, int> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            size.push_back(2);
            chunk.push_back(2);
            offset.push_back(0);
            ar.write(path, get_pointer(value), size, chunk, offset);
        }

        void load(
              archive & ar
            , std::string const & path
            , std::pair<int, int> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            chunk.push_back(2);
            offset.push_back(0);
            ar.read(path, get_pointer(value), chunk, offset);
        }

    }
}

int main() {

    std::string const filename = "example_pair_int_vectorizable.h5";

    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    std::pair<int, int> scalar_read, scalar_write(1, 2);
    std::vector<std::pair<int, int> > vector_read, vector_write(10, std::make_pair(3, 4));

    {
        alps::hdf5::archive ar(filename, "a");
        ar["/enum/scalar"] << scalar_write;
        ar["/enum/vector"] << vector_write;
    }

    {
        alps::hdf5::archive ar(filename);
        ar["/enum/scalar"] >> scalar_read;
        ar["/enum/vector"] >> vector_read;
    }

    boost::filesystem::remove(boost::filesystem::path(filename));

    bool match = (
           scalar_read == scalar_write
        && vector_read.size() == vector_write.size()
        && std::equal(vector_read.begin(), vector_read.end(), vector_write.begin())
    );
    std::cout << (match ? "SUCCESS" : "FAILURE") << std::endl;
    return match ? EXIT_SUCCESS : EXIT_FAILURE;
}
