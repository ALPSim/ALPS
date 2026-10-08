// SPDX-License-Identifier: MIT
#include <alps/ngs/params.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <complex>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template<class Action>
void rejects_storage(Action action, std::string const& path, bool scalar) {
    bool rejected = false;
    try { action(); }
    catch (std::runtime_error const& error) {
        std::string diagnostic(error.what());
        EXPECT_TRUE(diagnostic.find(scalar ? "Unsupported parameter scalar datatype"
                                       : "Unsupported parameter array datatype") != std::string::npos) << "Unsupported storage must have a specific diagnostic: " + diagnostic;
        EXPECT_TRUE(diagnostic.find(path) != std::string::npos) << "Unsupported storage diagnostic must identify its dataset: " + diagnostic;
        rejected = true;
    }
    EXPECT_TRUE(rejected) << "Unsupported parameter storage was accepted: " + path;
}

template<class T>
void unsupported_type(alps::hdf5::archive& archive, std::string const& name, T value) {
    std::string scalar_path = "/unsupported/" + name;
    archive[scalar_path] << value;
    archive.set_context(scalar_path);
    alps::detail::paramvalue existing(23);
    rejects_storage([&] { existing.load(archive); }, scalar_path, true);
    EXPECT_EQ(existing.cast<int>(), 23) << "Failed scalar load changed the existing value";

    std::string array_path = scalar_path + "_array";
    archive[array_path] << std::vector<T>{value, value};
    archive.set_context(array_path);
    rejects_storage([&] { existing.load(archive); }, array_path, false);
    EXPECT_EQ(existing.cast<int>(), 23) << "Failed array load changed the existing value";
}

template<class T>
void supported_value(alps::params& expected, alps::params const& actual,
                     std::string const& key, T const& value) {
    EXPECT_EQ(actual[key].cast<T>(), value) << "Supported checkpoint changed value: " + key;
    ASSERT_NE(actual.find(key), nullptr) << key;
    EXPECT_EQ(actual.find(key)->which(), expected.find(key)->which()) << "Supported checkpoint changed its native variant type: " + key;
}

void supported_checkpoint(alps::hdf5::archive& archive) {
    alps::params expected;
    expected["integer"] = -17;
    expected["real"] = 1.25;
    expected["boolean"] = true;
    expected["text"] = std::string("unchanged");
    expected["complex"] = std::complex<double>(2., -3.);
    expected["integers"] = std::vector<int>{-2, 0, 4};
    expected["reals"] = std::vector<double>{-1.5, 0., 2.25};
    expected["booleans"] = std::vector<bool>{true, false, true};
    expected["strings"] = std::vector<std::string>{"first", "", "last"};
    expected["complexes"] = std::vector<std::complex<double>>{{1., 2.}, {-3., 4.}};
    archive["/supported"] << expected;
    alps::params actual(archive, "/supported");
    EXPECT_EQ(actual.size(), expected.size()) << "Supported checkpoint lost parameters";
    supported_value(expected, actual, "integer", -17);
    supported_value(expected, actual, "real", 1.25);
    supported_value(expected, actual, "boolean", true);
    supported_value(expected, actual, "text", std::string("unchanged"));
    supported_value(expected, actual, "complex", std::complex<double>(2., -3.));
    supported_value(expected, actual, "integers", std::vector<int>{-2, 0, 4});
    supported_value(expected, actual, "reals", std::vector<double>{-1.5, 0., 2.25});
    supported_value(expected, actual, "booleans", std::vector<bool>{true, false, true});
    supported_value(expected, actual, "strings", std::vector<std::string>{"first", "", "last"});
    supported_value(expected, actual, "complexes", std::vector<std::complex<double>>{{1., 2.}, {-3., 4.}});
}

void transactional_reload(alps::hdf5::archive& archive) {
    archive["/transaction/a_supported"] << 10;
    archive["/transaction/z_unsupported"] << 1.25F;
    alps::params existing;
    existing["keep_first"] = std::string("original");
    existing["keep_second"] = 42;
    archive.set_context("/transaction");
    rejects_storage([&] { existing.load(archive); }, "/transaction/z_unsupported", true);
    EXPECT_TRUE(existing.size() == 2 && existing["keep_first"].cast<std::string>() == "original"
            && existing["keep_second"].cast<int>() == 42) << "Rejected checkpoint changed existing parameter values";
    ASSERT_EQ(existing.size(), 2u);
    auto key = existing.begin();
    EXPECT_EQ(key->first, "keep_first") << "Rejected checkpoint changed parameter iteration order";
    ++key;
    EXPECT_EQ(key->first, "keep_second") << "Rejected checkpoint changed parameter iteration order";
    EXPECT_TRUE(!existing.defined("a_supported") && !existing.defined("z_unsupported")) << "Rejected checkpoint published a partial result";

    rejects_storage([&] { alps::params rejected(archive, "/transaction"); },
                    "/transaction/z_unsupported", true);
}

void custom_reader(alps::hdf5::archive& archive) {
    archive["/custom/value"] << 1.25F;
    alps::params parameters;
    int calls = 0;
    parameters.set_value_reader([&](alps::hdf5::archive& reader) {
        ++calls;
        EXPECT_EQ(reader.get_context(), "/custom/value") << "Custom reader receives the dataset context";
        float value = 0;
        reader[""] >> value;
        return alps::detail::paramvalue(static_cast<double>(value));
    });
    archive.set_context("/custom");
    parameters.load(archive);
    EXPECT_TRUE(calls == 1 && parameters["value"].cast<double>() == 1.25) << "Custom reader must retain responsibility for its own decoding";
    EXPECT_EQ(archive.get_context(), "/custom") << "Successful custom load changed the caller context";
    parameters.load(archive);
    EXPECT_TRUE(calls == 2 && parameters["value"].cast<double>() == 1.25) << "Reload must preserve the configured custom reader";
}
}

class ParamsCheckpoint : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive archive{(temporary.path() / "checkpoint.h5").string(), "w"};
};
TEST_F(ParamsCheckpoint, SupportedValuesPreserveTypes) { supported_checkpoint(archive); }
TEST_F(ParamsCheckpoint, RejectsUnsupportedFloatStorage) { unsupported_type(archive, "float", 1.25F); }
TEST_F(ParamsCheckpoint, RejectsUnsupportedUnsignedStorage) { unsupported_type(archive, "unsigned", 42U); }
TEST_F(ParamsCheckpoint, RejectsUnsupportedWideIntegerStorage) {
    if (sizeof(long long) <= sizeof(int)) GTEST_SKIP() << "long long uses native int storage";
    unsupported_type(archive, "wide_integer", 1099511627776LL);
}
TEST_F(ParamsCheckpoint, RejectsUnsupportedLongStorage) {
    // On LLP64, native long has int's storage and remains supported.
    if (sizeof(long) <= sizeof(int)) GTEST_SKIP() << "long uses native int storage";
    unsupported_type(archive, "native_long", static_cast<long>(1099511627776LL));
}
TEST_F(ParamsCheckpoint, FailedReloadIsTransactional) { transactional_reload(archive); }
TEST_F(ParamsCheckpoint, CustomReaderSurvivesReload) { custom_reader(archive); }
