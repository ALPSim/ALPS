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

typedef enum { PLUS, MINUS } enum_type;

namespace alps {
    namespace hdf5 {

        template<> struct scalar_type<enum_type> {
            typedef boost::int_t<sizeof(enum_type) * 8>::exact type;
        };

        template<> struct is_continuous<enum_type> : public boost::true_type {};

        template<> struct has_complex_elements<enum_type> : public boost::false_type {};

        namespace detail {

            template<> struct get_extent<enum_type> {
                static std::vector<std::size_t> apply(enum_type const & value) {
                    return std::vector<std::size_t>();
                }
            };

            template<> struct set_extent<enum_type> {
                static void apply(enum_type &, std::vector<std::size_t> const &) {}
            };

            template<> struct is_vectorizable<enum_type> {
                static bool apply(enum_type const & value) {
                    return true;
                }
            };

            template<> struct get_pointer<enum_type> {
                static scalar_type<enum_type>::type * apply(enum_type & value) {
                    return reinterpret_cast<scalar_type<enum_type>::type *>(&value);
                }
            };
        
            template<> struct get_pointer<enum_type const> {
                static scalar_type<enum_type>::type const * apply(enum_type const & value) {
                    return reinterpret_cast<scalar_type<enum_type>::type const *>(&value);
                }
            };
        }

        void save(
              alps::hdf5::archive & ar
            , std::string const & path
            , enum_type const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (size.size() == 0) {
                size.push_back(1);
                chunk.push_back(1);
                offset.push_back(0);
            }
            ar.write(path, (scalar_type<enum_type>::type const *)get_pointer(value), size, chunk, offset);
        }
        
        void load(
              alps::hdf5::archive & ar
            , std::string const & path
            , enum_type & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (chunk.size() == 0) {
                chunk.push_back(1);
                offset.push_back(0);
            }
            ar.read(path, (scalar_type<enum_type>::type *)get_pointer(value), chunk, offset);
        }

    }
}

int main() {

    std::string const filename = "example_enum_vectorizable.h5";

    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    enum_type scalar_read, scalar_write = PLUS;
    std::vector<enum_type> vector_read, vector_write(10, MINUS);

    {
        alps::hdf5::archive ar(filename, "a");
        ar 
          << alps::make_pvp("/enum/scalar", scalar_write)
          << alps::make_pvp("/enum/vector", vector_write)
        ;
    }

    {
        alps::hdf5::archive ar(filename);
        ar 
          >> alps::make_pvp("/enum/scalar", scalar_read)
          >> alps::make_pvp("/enum/vector", vector_read)
        ;
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
