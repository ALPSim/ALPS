#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/ngs/params.hpp>
#include <sstream>

namespace {
template<class Parameters> void expect_missing(Parameters& parameters) {
    std::ostringstream output;
    try {
        output << parameters["not_in_parms"];
        FAIL() << "Missing parameter did not throw";
    } catch (const std::exception& error) {
        const std::string diagnostic(error.what());
        EXPECT_EQ(diagnostic.substr(0, diagnostic.find('\n')), "No parameter 'not_in_parms' available");
    }
}
}
TEST(ParamsLookup, MissingMutableAndConstValuesReportTheirKey) {
    alps::params parameters;
    parameters["hello"] = "world";
    EXPECT_EQ(parameters["hello"].cast<std::string>(), "world");
    expect_missing(parameters);
    const alps::params copy(parameters);
    expect_missing(copy);
}
