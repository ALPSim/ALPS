/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/scheduler/montecarlo.h>
#include <alps/parser/xmlparser.h>
#include <alps/parser/xmlstream.h>
#include <gtest/gtest.h>
#include <sstream>
#include <valarray>

namespace {

// Results written as XML and read back, by the reader of the task files of
// the scheduler (`read_xml_averages`) and that of parapack (the handler).
std::vector<alps::alea::result_set> round_trip(const alps::alea::result_set& results)
{
  std::ostringstream text;
  {
    alps::oxstream out(text);
    alps::scheduler::write_xml_averages(out, results);
  }
  std::vector<alps::alea::result_set> sets;
  alps::scheduler::AveragesXMLHandler handler(sets);
  alps::XMLParser parser(handler);
  std::istringstream sax(text.str());
  parser.parse(sax);
  std::istringstream stream(text.str());
  const alps::XMLTag tag = alps::parse_tag(stream);
  sets.push_back(alps::scheduler::read_xml_averages(stream, tag));
  return sets;
}

// XML keeps the error to three digits and the mean to four digits more than
// the error
void expect_close(const alps::alea::observable_result& written, const alps::alea::observable_result& read)
{
  ASSERT_EQ(read.size(), written.size());
  EXPECT_EQ(read.count(), written.count());
  EXPECT_EQ(read.is_scalar(), written.is_scalar());
  EXPECT_EQ(read.sign(), written.sign());
  EXPECT_EQ(read.labels(), written.labels());
  EXPECT_FALSE(read.has_bins());
  for (std::size_t i = 0; i != written.size(); ++i) {
    SCOPED_TRACE(i);
    const double error = written.error_vector()[i];
    EXPECT_NEAR(read.mean_vector()[i], written.mean_vector()[i], 1e-3 * error);
    EXPECT_NEAR(read.error_vector()[i], error, 5e-3 * error);
    EXPECT_EQ(read.converged_vector()[i], written.converged_vector()[i]);
    const std::vector<double> levels = written.binning_errors(i), read_levels = read.binning_errors(i);
    ASSERT_EQ(read_levels.size(), levels.size());
    for (std::size_t level = 0; level != levels.size(); ++level)
      EXPECT_NEAR(read_levels[level], levels[level], 5e-3 * levels[level]);
    const double tau = written.tau_vector()[i], variance = written.variance_vector()[i];
    if (tau == tau)
      EXPECT_NEAR(read.tau_vector()[i], tau, 5e-3 * std::abs(tau));
    if (variance == variance)
      EXPECT_NEAR(read.variance_vector()[i], variance, 5e-3 * variance);
  }
}

} // namespace

TEST(XmlAverages, ResultsComeBackAsWritten) {
  alps::alea::observable_set set;
  set << alps::alea::observable("Energy") << alps::alea::observable("Sign")
      << alps::alea::observable("Density").signed_by("Sign")
      << alps::alea::observable("Correlations").labels({"( 0 )", "( 1 )"})
      << alps::alea::observable("Bond-type Energy")
      << alps::alea::observable("Counter").without_bins()
      << alps::alea::histogram("Steps", -2, 2);
  double x = 0.3;
  for (int i = 0; i != 4096; ++i) {
    x = 0.9 * x + 0.1 * std::sin(1.7 * i);     // correlated, reproducible
    set["Energy"] << -2. + x;
    set["Sign"] << (i % 7 == 0 ? -1. : 1.);
    set["Density"] << (i % 7 == 0 ? -1. : 1.) * (0.5 + 0.1 * x);
    set["Correlations"] << std::vector<double>{1. + x, x * x};
    set["Bond-type Energy"] << std::valarray<double>(-2. + x, 1);
    set["Counter"] << double(i % 3);
    set["Steps"] << double(i % 5 - 2);
  }
  alps::alea::result_set results = set.results();
  results.insert("Energy^2", results["Energy"] * results["Energy"]);

  for (const alps::alea::result_set& read : round_trip(results)) {
    ASSERT_EQ(read.names(), results.names());
    for (const std::string& name : results.names()) {
      SCOPED_TRACE(name);
      if (name == "Steps") {
        EXPECT_TRUE(read[name].is_histogram());
        EXPECT_EQ(read[name].histogram(), results[name].histogram());
      } else {
        expect_close(results[name], read[name]);
      }
    }
  }
}

// As ALPS 3.0 wrote them, with the means of every binning level
TEST(XmlAverages, ReadsTheAveragesOfAlps30) {
  std::istringstream text(R"(<AVERAGES>
    <SCALAR_AVERAGE name="Energy">
      <COUNT>20000</COUNT>
      <MEAN method="simple">-36.4565</MEAN>
      <ERROR converged="yes" method="binning">0.247</ERROR>
      <VARIANCE method="simple">80.7</VARIANCE>
      <AUTOCORR method="binning">7.08</AUTOCORR>
      <BINNED size="1"><COUNT>20000</COUNT><MEAN method="simple">-36.4565</MEAN><ERROR method="simple">0.0635</ERROR></BINNED>
      <BINNED size="2"><COUNT>10000</COUNT><MEAN method="simple">-36.4565</MEAN><ERROR method="simple">0.0873</ERROR></BINNED>
    </SCALAR_AVERAGE>
    <VECTOR_AVERAGE name="Bond-type Energy" nvalues="1">
      <SCALAR_AVERAGE indexvalue="0">
        <COUNT>20000</COUNT>
        <MEAN method="simple">-36.4565</MEAN>
        <ERROR converged="maybe" method="binning">0.247</ERROR>
      </SCALAR_AVERAGE>
    </VECTOR_AVERAGE>
    <SCALAR_AVERAGE name="Specific Heat" signed="true">
      <COUNT>20000</COUNT>
      <MEAN method="jackknife">0.77446</MEAN>
      <ERROR converged="no" method="jackknife">0.0113</ERROR>
    </SCALAR_AVERAGE>
  </AVERAGES>)");
  const alps::XMLTag tag = alps::parse_tag(text);
  const alps::alea::result_set read = alps::scheduler::read_xml_averages(text, tag);
  ASSERT_EQ(read.names(), (std::vector<std::string>{"Bond-type Energy", "Energy", "Specific Heat"}));

  const alps::alea::observable_result& energy = read["Energy"];
  EXPECT_TRUE(energy.is_scalar());
  EXPECT_EQ(energy.count(), 20000u);
  EXPECT_DOUBLE_EQ(energy.mean(), -36.4565);
  EXPECT_DOUBLE_EQ(energy.error(), 0.247);
  EXPECT_DOUBLE_EQ(energy.tau(), 7.08);
  EXPECT_DOUBLE_EQ(energy.variance_vector()[0], 80.7);
  EXPECT_EQ(energy.converged(), alps::alea::CONVERGED);
  EXPECT_EQ(energy.binning_errors(), (std::vector<double>{0.0635, 0.0873}));

  const alps::alea::observable_result& bond = read["Bond-type Energy"];
  EXPECT_FALSE(bond.is_scalar());
  EXPECT_EQ(bond.size(), 1u);
  EXPECT_TRUE(bond.labels().empty());    // an index is no label
  EXPECT_EQ(bond.converged_vector()[0], alps::alea::MAYBE_CONVERGED);
  EXPECT_TRUE(std::isnan(bond.tau_vector()[0]));

  EXPECT_EQ(read["Specific Heat"].sign(), "Sign");
  EXPECT_EQ(read["Specific Heat"].converged(), alps::alea::NOT_CONVERGED);
}
