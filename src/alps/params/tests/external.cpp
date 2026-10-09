// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/ngs/params.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <memory>
#include <sstream>
#include <stdexcept>

struct source final : alps::detail::paramvalue_source {
    std::vector<double> values{1., 2.};
    alps::detail::paramvalue native_value() const override { return values; }
    void save(alps::hdf5::archive & ar) const override { ar[""] << values; }
    void print(std::ostream & out) const override { out << values.front(); }
    void * object(char const *) const override { return nullptr; }
};

TEST(ParamsExternal, SharedLifetimeAndNativeSerialization) {
    auto value = std::make_shared<source>();
    std::weak_ptr<source> lifetime = value;
    alps::params parameters;
    parameters["vector"] = alps::detail::paramvalue(value);
    alps::params copy(parameters);
    value->values[0] = 9.;
    EXPECT_EQ(copy["vector"].cast<std::vector<double>>()[0], 9.);
    EXPECT_EQ(parameters.find("vector")->cast<std::vector<double>>()[0], 9.);

    // Boost serialization materializes native values instead of attempting
    // to serialize a language runtime pointer or callback.
    std::stringstream buffer;
    { boost::archive::text_oarchive archive(buffer); archive << parameters; }
    alps::params restored;
    { boost::archive::text_iarchive archive(buffer); archive >> restored; }
    EXPECT_EQ(restored["vector"].cast<std::vector<double>>(), value->values);

    parameters["vector"] = 3;
    EXPECT_EQ(parameters["vector"].cast<int>(), 3);
    value.reset();
    EXPECT_TRUE(!lifetime.expired());
    copy.erase("vector");
    EXPECT_TRUE(lifetime.expired());

}

TEST(ParamsConversion, IntegerOverflowAndStringLists) {
    alps::params parameters;
    parameters["wide"] = std::string("9007199254740993");
    EXPECT_EQ(parameters["wide"].cast<long long>(), 9007199254740993LL);
    EXPECT_THROW(parameters["wide"].cast<int>(), std::out_of_range);

    parameters["names"] = std::vector<std::string>{"Energy", "Stiffness"};
    EXPECT_EQ(parameters["names"].cast<std::string>(), "Energy,Stiffness");
    parameters["names"] = std::vector<std::string>{"", "middle", ""};
    EXPECT_EQ(parameters["names"].cast<std::string>(), ",middle,");

}

TEST(ParamsExternal, HeterogeneousArchiveListMaterializesNativeValues) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive archive((temporary.path() / "external_list.h5").string(), "w");
    archive["/list/0"] << true;
    archive["/list/1"] << 2;
    archive["/list/2"] << 10.5;
    archive.set_context("/list");
    alps::detail::paramvalue list;
    list.load(archive);
    EXPECT_EQ(list.cast<std::vector<double>>(), std::vector<double>({1., 2., 10.5}));
    EXPECT_EQ(list.cast<std::vector<int>>(), std::vector<int>({1, 2, 10}));
    list.save(archive);
    list.load(archive);
    EXPECT_EQ(list.cast<std::vector<int>>(), std::vector<int>({1, 2, 10}));
}
