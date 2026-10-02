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

        paramvalue_source::~paramvalue_source() = default;

        namespace {
            // Retain scalar types when a Python list checkpoint is resumed by
            // a native simulation. Reuse the existing value-provider contract
            // so conversion happens in the type requested by the consumer.
            class checkpoint_list final : public paramvalue_source {
            public:
                explicit checkpoint_list(std::vector<paramvalue> values)
                    : values_(std::move(values)) {}
                bool native_elements(std::vector<paramvalue> & values) const override {
                    values = values_;
                    return true;
                }
                paramvalue native_value() const override {
                    bool text = false, complex = false, real = false, integer = false;
                    for (auto const & value : values_) {
                        text |= value.which() == paramvalue_index<std::string>::value;
                        complex |= value.which() == paramvalue_index<std::complex<double>>::value;
                        real |= value.which() == paramvalue_index<double>::value;
                        integer |= value.which() == paramvalue_index<int>::value;
                    }
                    if (text) return converted<std::string>();
                    if (complex) return converted<std::complex<double>>();
                    if (real) return converted<double>();
                    if (integer) return converted<int>();
                    return converted<bool>();
                }
                void save(hdf5::archive & ar) const override {
                    if (ar.is_data("")) ar.delete_data("");
                    if (ar.is_group("")) ar.delete_group("");
                    ar.create_group("");
                    for (std::size_t i = 0; i < values_.size(); ++i)
                        ar[std::to_string(i)] << values_[i];
                }
                void print(std::ostream & out) const override { out << native_value(); }
                void * object(char const *) const override { return nullptr; }
            private:
                template <typename T> std::vector<T> converted() const {
                    std::vector<T> values;
                    values.reserve(values_.size());
                    for (auto const & value : values_) values.push_back(value.cast<T>());
                    return values;
                }
                std::vector<paramvalue> values_;
            };
        }

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
                return cast<T>();                                                   \
            }
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_OPERATOR_T_IMPL)
        #undef ALPS_NGS_PARAMVALUE_OPERATOR_T_IMPL

        #define ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL(T)                                \
            paramvalue & paramvalue::operator=( T const & arg) {                    \
                paramvalue_base::operator=(arg);                                    \
                source_.reset();                                                   \
                return *this;                                                        \
            }
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL)
        #undef ALPS_NGS_PARAMVALUE_OPERATOR_EQ_IMPL

        void paramvalue::save(hdf5::archive & ar) const {
            if (source_) {
                source_->save(ar);
                return;
            }
            boost::apply_visitor(
                paramvalue_saver(ar), static_cast<paramvalue_base const &>(*this)
            );
        }

        void paramvalue::load(hdf5::archive & ar) {
            if (ar.is_group("")) {
                // The Python archive stores Boolean/mixed lists as numbered
                // scalar children. Native-only simulations must be able to
                // resume those parameters too, without a Python decoder.
                auto children = ar.list_children("");
                if (children.empty())
                    throw std::runtime_error("Cannot load an empty parameter group" + ALPS_STACKTRACE);
                std::vector<paramvalue> values;
                values.reserve(children.size());
                for (std::size_t i = 0; i < children.size(); ++i) {
                    std::string child = std::to_string(i);
                    if (!ar.is_data(child))
                        throw std::runtime_error("Parameter group is not a scalar list" + ALPS_STACKTRACE);
                    if (ar.is_complex(child) && ar.dimensions(child) < 2) {
                        std::complex<double> number;
                        ar[child] >> number;
                        values.emplace_back(number);
                    } else if (!ar.is_scalar(child)) {
                        throw std::runtime_error("Parameter list contains a nonscalar value" + ALPS_STACKTRACE);
                    } else if (ar.is_datatype<signed char>(child)) {
                        // bool and int8 share their HDF5 storage type. Reading
                        // the byte directly preserves both 0/1 and signed data.
                        signed char number;
                        ar[child] >> number;
                        std::string type;
                        if (ar.is_attribute(child + "/@__alps_type__"))
                            ar[child + "/@__alps_type__"] >> type;
                        if (type == "int8") values.emplace_back(static_cast<int>(number));
                        else values.emplace_back(number != 0);
                    } else if (ar.is_datatype<double>(child)
                               || ar.is_datatype<float>(child)
                               || ar.is_datatype<long double>(child)) {
                        double number;
                        ar[child] >> number;
                        values.emplace_back(number);
                    } else if (ar.is_datatype<int>(child)) {
                        int number;
                        ar[child] >> number;
                        values.emplace_back(number);
                    } else {
                        // Keep wider integers exact despite the native
                        // variant's int-sized integer alternative.
                        std::string value;
                        ar[child] >> value;
                        values.emplace_back(value);
                    }
                }
                *this = paramvalue(std::make_shared<checkpoint_list>(std::move(values)));
                return;
            }
            #define ALPS_NGS_PARAMVALUE_LOAD_HDF5(T)                                \
                {                                                                    \
                    T value;                                                        \
                    ar[""] >> value;                                        \
                    operator=(value);                                                \
                }
            #define ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(T, U)                        \
                else if (ar.is_datatype< T >(""))                                    \
                    ALPS_NGS_PARAMVALUE_LOAD_HDF5(U)
            // A complex scalar is stored as a trailing dimension of two
            // reals, so archive::is_scalar() reports false for it and it fell
            // into the vector branch below, where loading it as
            // vector<complex> failed with "dimensions do not match". Rank
            // tells them apart: a complex scalar has dimensions() == 1, a
            // complex vector -- even a one-element one -- has 2.
            if (ar.is_complex("") && ar.dimensions("") < 2)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5(std::complex<double>)
            else if (ar.is_scalar("")) {
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
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(bool, std::vector<bool>)
                ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK(
                    std::string, std::vector<std::string>
                )
            }
            #undef ALPS_NGS_PARAMVALUE_LOAD_HDF5
            #undef ALPS_NGS_PARAMVALUE_LOAD_HDF5_CHECK
        }

        std::ostream & operator<<(std::ostream & os, paramvalue const & arg) {
            if (arg.source()) {
                arg.source()->print(os);
                return os;
            }
            paramvalue_ostream visitor(os);
            boost::apply_visitor(visitor, arg);
            return os;
        }        
    }
}
