/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_PARAMVALUE_HPP
#define ALPS_NGS_DETAIL_PARAMVALUE_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/detail/remove_cvr.hpp>
#include <alps/ngs/detail/paramvalue_reader.hpp>

#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp> 
#include <boost/serialization/complex.hpp> 
#include <boost/serialization/split_member.hpp>

#include <string>
#include <complex>
#include <ostream>
#include <stdexcept>
#include <memory>
#include <utility>

#define ALPS_NGS_FOREACH_PARAMETERVALUE_ADDABLE_TYPE(CALLBACK)                      \
    CALLBACK(double)                                                                \
    CALLBACK(int)                                                                   \
    CALLBACK(bool)                                                                  \
    CALLBACK(std::string)                                                           \
    CALLBACK(std::complex<double>)                                                  \
    CALLBACK(std::vector<double>)                                                   \
    CALLBACK(std::vector<int>)                                                      \
    CALLBACK(std::vector<std::string>)                                              \
    CALLBACK(std::vector<std::complex<double> >)

#define ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(CALLBACK)                    \
    ALPS_NGS_FOREACH_PARAMETERVALUE_ADDABLE_TYPE(CALLBACK)                          \
    CALLBACK(std::vector<bool>)

    #define ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(CALLBACK)                          \
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(CALLBACK)

namespace alps {

    namespace detail {

		template <typename T> struct paramvalue_index {};
		template <> struct paramvalue_index<double> {
			enum { value = 0 };
		};
		template <> struct paramvalue_index<int> {
			enum { value = 1 };
		};
		template <> struct paramvalue_index<bool> {
			enum { value = 2 };
		};
		template <> struct paramvalue_index<std::string> {
			enum { value = 3 };
		};
		template <> struct paramvalue_index<std::complex<double> > {
			enum { value = 4 };
		};
		template <> struct paramvalue_index<std::vector<double> > {
			enum { value = 5 };
		};
		template <> struct paramvalue_index<std::vector<int> > {
			enum { value = 6 };
		};
		template <> struct paramvalue_index<std::vector<std::string> > {
			enum { value = 7 };
		};
		template <> struct paramvalue_index<std::vector<std::complex<double> > > {
			enum { value = 8 };
		};
		template <> struct paramvalue_index<std::vector<bool> > {
			enum { value = 9 };
		};

        class paramvalue;

        // Optional value provider for language bindings. No interpreter
        // headers or runtime are required by the native library. A binding
        // can retain a mutable value and supply a checked native snapshot
        // only when a C++ consumer actually requests it.
        struct ALPS_DECL paramvalue_source {
            virtual ~paramvalue_source();
            virtual paramvalue native_value() const = 0;
            // A heterogeneous sequence must be converted element by element
            // to the type requested by the consumer, without first forcing
            // every element into one variant alternative.
            virtual bool native_elements(std::vector<paramvalue> &) const { return false; }
            virtual void save(hdf5::archive &) const = 0;
            virtual void print(std::ostream &) const = 0;
            virtual void * object(char const * binding) const = 0;
            // Some bindings give sequences a specific textual form. Keep
            // that conversion separate from a homogeneous numeric snapshot.
            virtual bool native_text(std::string &) const { return false; }
        };

    }

    template<typename T> T extract (detail::paramvalue const & arg);

    namespace detail {

        template <typename T> struct paramvalue_vector : std::false_type {};
        template <typename T> struct paramvalue_vector<std::vector<T>> : std::true_type {};

        template<class Archive> struct paramvalue_serializer 
            : public boost::static_visitor<> 
        {
            public:

                paramvalue_serializer(Archive & a)
                    : ar(a)
                {}

                template <typename U> void operator()(U & v) const {
                    std::size_t type = paramvalue_index<typename remove_cvr<U>::type>::value;
                    ar
                        << type
                        << v
                    ;
                }

            private:

                Archive & ar;
        };

        #define ALPS_NGS_PARAMVALUE_VARIANT_TYPE(T)    T,
        typedef boost::mpl::pop_back<boost::mpl::vector<
            ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(
                ALPS_NGS_PARAMVALUE_VARIANT_TYPE
            ) void
        >::type >::type paramvalue_types;
        #undef ALPS_NGS_PARAMVALUE_VARIANT_TYPE
        typedef boost::make_variant_over<paramvalue_types>::type paramvalue_base;

        class paramvalue : public paramvalue_base {

            public:

                paramvalue() {}

                paramvalue(paramvalue const & v)
                    : paramvalue_base(static_cast<paramvalue_base const &>(v))
                    , source_(v.source_)
                {}

                explicit paramvalue(std::shared_ptr<paramvalue_source> source)
                    : source_(std::move(source)) {}

                std::shared_ptr<paramvalue_source> const & source() const { return source_; }

                paramvalue const& operator=(paramvalue const& x)
                {
                  static_cast<paramvalue_base&>(*this) = static_cast<paramvalue_base const&>(x);
                  source_ = x.source_;
                  return *this;
                }
                template<typename T> T cast() const {
                    if (source_) {
                        if constexpr (std::is_same<T, std::string>::value) {
                            std::string value;
                            if (source_->native_text(value))
                                return value;
                        }
                        if constexpr (paramvalue_vector<T>::value) {
                            std::vector<paramvalue> elements;
                            if (source_->native_elements(elements)) {
                                T values;
                                values.reserve(elements.size());
                                for (auto const & element : elements)
                                    values.push_back(element.cast<typename T::value_type>());
                                return values;
                            }
                        }
                        return source_->native_value().cast<T>();
                    }
                    paramvalue_reader< T > visitor;
                    boost::apply_visitor(visitor, *this);
                    return visitor.get_value();
                }

                #define ALPS_NGS_PARAMVALUE_MEMBER_DECL(T)                          \
                    paramvalue( T const & v) : paramvalue_base(v) {}                \
                    operator T () const;                                            \
                    paramvalue & operator=( T const &);
                ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_MEMBER_DECL)
                #undef ALPS_NGS_PARAMVALUE_MEMBER_DECL

                
                #define ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(T)                  \
                    paramvalue(T const & value)                                     \
                        : paramvalue_base(static_cast<int>(value))                  \
                    {}                                                              \
                    paramvalue & operator=(T const & value) {                       \
                        return *this = static_cast<int>(value);                     \
                    }
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(char)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(unsigned char)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(short)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(unsigned short)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(unsigned)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(long)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(unsigned long)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(long long)
                ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT(unsigned long long)
                #undef ALPS_NGS_PARAMVALUE_MEMBER_DECL_CONVERT

                paramvalue(long double const & value)
                    : paramvalue_base(static_cast<double>(value))
                {}
                paramvalue & operator=(long double const & value) {
                    return *this = static_cast<double>(value);
                }

                paramvalue(char const * v) : paramvalue_base(std::string(v)) {}
                paramvalue & operator=(char const * value) { 
                    return *this = std::string(value); 
                }

                void save(hdf5::archive &) const;
                void load(hdf5::archive &);
                
            private:

                std::shared_ptr<paramvalue_source> source_;
            
                friend class boost::serialization::access;

                template<class Archive> void save(
                    Archive & ar, const unsigned int version
                ) const {
                    paramvalue_serializer<Archive> visitor(ar);
                    if (source_) {
                        paramvalue native = source_->native_value();
                        boost::apply_visitor(visitor, native);
                    } else
                        boost::apply_visitor(visitor, *this);
                }

                template<class Archive> void load(
                    Archive & ar, const unsigned int
                ) {
                    std::size_t type;
                    ar >> type;
                    if (false);
                    #define ALPS_NGS_PARAMVALUE_LOAD(T)                                 \
                        else if (type == paramvalue_index< T >::value) {                \
                            T value;                                                    \
                            ar >> value;                                                \
                            operator=(value);                                           \
                        }
                    // TODO: fix python serialization!
                    // ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMVALUE_LOAD)
                    ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE_NO_PYTHON(ALPS_NGS_PARAMVALUE_LOAD)
                    #undef ALPS_NGS_PARAMVALUE_LOAD
                    else
                        throw std::runtime_error("unknown type" + ALPS_STACKTRACE);
                }

                BOOST_SERIALIZATION_SPLIT_MEMBER()
        };

        ALPS_DECL std::ostream & operator<<(std::ostream & os, paramvalue const & arg);

        template<typename T> T extract_impl (paramvalue const & arg, T) {
            return arg.cast<T>();
        }
    }

    template<typename T> struct cast_hook<T, detail::paramvalue> {
        static inline std::complex<T> apply(detail::paramvalue const & arg) {
            return extract<T>(arg);
        }
    };

    template<typename T> T extract (detail::paramvalue const & arg) {
        return arg.cast<T>();
    }
}

#endif
