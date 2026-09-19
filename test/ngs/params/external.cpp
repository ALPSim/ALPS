// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/ngs/params.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <memory>
#include <sstream>
#include <stdexcept>

void require(bool value) {
    if (!value) throw std::runtime_error("external parameter contract failed");
}

struct source final : alps::detail::paramvalue_source {
    std::vector<double> values{1., 2.};
    alps::detail::paramvalue native_value() const override { return values; }
    void save(alps::hdf5::archive & ar) const override { ar[""] << values; }
    void print(std::ostream & out) const override { out << values.front(); }
    void * object(char const *) const override { return nullptr; }
};

int main() {
    auto value = std::make_shared<source>();
    std::weak_ptr<source> lifetime = value;
    alps::params parameters;
    parameters["vector"] = alps::detail::paramvalue(value);
    alps::params copy(parameters);
    value->values[0] = 9.;
    require(copy["vector"].cast<std::vector<double>>()[0] == 9.);
    require(parameters.find("vector")->cast<std::vector<double>>()[0] == 9.);

    // Boost serialization materializes native values instead of attempting
    // to serialize a language runtime pointer or callback.
    std::stringstream buffer;
    { boost::archive::text_oarchive archive(buffer); archive << parameters; }
    alps::params restored;
    { boost::archive::text_iarchive archive(buffer); archive >> restored; }
    require(restored["vector"].cast<std::vector<double>>() == value->values);

    parameters["vector"] = 3;
    require(parameters["vector"].cast<int>() == 3);
    value.reset();
    require(!lifetime.expired());
    copy.erase("vector");
    require(lifetime.expired());

    parameters["wide"] = std::string("9007199254740993");
    require(parameters["wide"].cast<long long>() == 9007199254740993LL);
    try {
        parameters["wide"].cast<int>();
        throw std::runtime_error("overflowing conversion unexpectedly succeeded");
    } catch (std::out_of_range const &) {}
}
