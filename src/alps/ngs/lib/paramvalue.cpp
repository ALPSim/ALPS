/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/short_print.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/ngs/detail/paramvalue.hpp>
#include <alps/ngs/detail/type_wrapper.hpp>

namespace alps {
    namespace detail {

        struct paramvalue_saver: public boost::static_visitor<> {

            paramvalue_saver(hdf5::archive & a)
                : ar(a) 
            {}

            template<typename T> void operator()(T const & v) const {
                ar[""] << v;
            }

            hdf5::archive & ar;
        };

        struct paramvalue_ostream : public boost::static_visitor<> {
            public:

                paramvalue_ostream(std::ostream & arg) : os(arg) {}

                template <typename U> void operator()(U const & v) const {
                    os << short_print(v);
                }

            private:

                std::ostream & os;
        };

        #define ALPS_NGS_PARAMVALUE_OPERATOR_T_IMPL(T)                                \
            paramvalue::operator T () const {                                        \
                paramvalue_reader< T > visitor;                                        \
                boost::apply_visitor(visitor, *this);                               \
                return visitor.get_value();                                            \
            }
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_OPERATOR_T_IMPL)
        #undef ALPS_NGS_PARAMVALUE_OPERATOR_T_IMPL

        #define ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL(T)                                \
            paramvalue & paramvalue::operator=( T const & arg) {                    \
                paramvalue_base::operator=(arg);                                    \
                return *this;                                                        \
            }
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL)
        #undef ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL

        void paramvalue::save(hdf5::archive & ar) const {
            boost::apply_visitor(
                paramvalue_saver(ar), static_cast<paramvalue_base const &>(*this)
            );
        }

        void paramvalue::load(hdf5::archive & ar) {
            #define ALPS_NGS_PARAMVALUE_LOAD_HDF5(T)                                \
                {                                                                    \
                    T value;                                                        \
                    ar[""] >> value;                                        \
                    operator=(value);                                                \
                }
            #define ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(T, U)                        \
                else if (ar.is_datatype< T >(""))                                    \
                    ALPS_NGS_PARAMVALUE_LOAD_HDF5(U)
            if (ar.is_scalar("")) {
                if (ar.is_complex(""))
                    ALPS_NGS_PARAMVALUE_LOAD_HDF5(std::complex<double>)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(double, double)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(int, int)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(bool, bool)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(std::string, std::string)
            } else {
                if (ar.is_complex(""))
                    ALPS_NGS_PARAMVALUE_LOAD_HDF5(
                        std::vector<std::complex<double> >
                    )
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(double, std::vector<double>)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(int, std::vector<int>)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(
                    std::string, std::vector<std::string>
                )
            }
            #undef ALPS_NGS_PARAMVALUE_LOAD_HDF5
            #undef ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK
        }

        std::ostream & operator<<(std::ostream & os, paramvalue const & arg) {
            paramvalue_ostream visitor(os);
            boost::apply_visitor(visitor, arg);
            return os;
        }        
    }
}
