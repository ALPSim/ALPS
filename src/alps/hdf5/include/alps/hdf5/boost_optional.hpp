/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Copied from ALPSCore 7146b9e1 hdf5/include/alps/hdf5/boost_optional.hpp
// (added in ALPSCore 7b5b7be5). Changes: it includes <alps/hdf5/archive.hpp>
// and the headers it uses, load() reads into a temporary, so loading into
// an empty optional is defined, and both functions reject non-native payload
// types with a static_assert. The empty flag is stored as the attribute
// "<path>/@alps_hdf5_optional_empty", as in ALPSCore, so files interchange.


#ifndef ALPS_HDF5_BOOST_OPTIONAL_4ea49b8c8c814a52b0b3c6ac5bcc7839
#define ALPS_HDF5_BOOST_OPTIONAL_4ea49b8c8c814a52b0b3c6ac5bcc7839

#include <alps/hdf5/archive.hpp>

#include <boost/optional.hpp>

#include <string>
#include <vector>

namespace alps {
    namespace hdf5 {

        template <typename T>
        void save(alps::hdf5::archive& ar, const std::string& path,
                  const boost::optional<T>& value,
                  std::vector<std::size_t> /*size*/=std::vector<std::size_t>(),
                  std::vector<std::size_t> /*chunk*/=std::vector<std::size_t>(),
                  std::vector<std::size_t> /*offset*/=std::vector<std::size_t>())
        {
            static_assert(is_native_type<T>::value,
                "boost::optional<T> is stored as one HDF5 value: T must be a native type");
            if (ar.is_group(path)) ar.delete_group(path);
            if (!value) {
                ar.write(path,T()); // FIXME?: T must have a default constructor then.
                ar.write(path+"/@alps_hdf5_optional_empty", true);
            } else {
                ar.write(path,*value);
                ar.write(path+"/@alps_hdf5_optional_empty", false);
            }
        }
          
        template <typename T>
        void load(alps::hdf5::archive& ar, const std::string& path,
                  boost::optional<T>& value,
                  std::vector<std::size_t> /*size*/=std::vector<std::size_t>(),
                  std::vector<std::size_t> /*chunk*/=std::vector<std::size_t>(),
                  std::vector<std::size_t> /*offset*/=std::vector<std::size_t>())
        {
            static_assert(is_native_type<T>::value,
                "boost::optional<T> is stored as one HDF5 value: T must be a native type");
            bool is_empty=false;
            ar.read(path+"/@alps_hdf5_optional_empty", is_empty);
            if (is_empty) {
                value = boost::none;
            } else {
                // ALPSCore reads into *value, which is undefined for an empty optional.
                T stored;
                ar.read(path, stored);
                value = stored;
            }
        }

    } // hdf5::
} // alps::

#endif /* ALPS_HDF5_BOOST_OPTIONAL_4ea49b8c8c814a52b0b3c6ac5bcc7839 */
