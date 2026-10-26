/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Regressions for the legacy alea binning, histogram and mcdata containers.

#include <gtest/gtest.h>
#include <alps/alea.h>
#include <alps/alea/histogramdata.h>
#include <alps/alea/mcdata.hpp>
#include <alps/parser/parser.h>
#include <sstream>
#include <stdexcept>

namespace {

alps::RealObservable one_to_eight()
{
  alps::RealObservable obs("x");
  for (int i = 1; i <= 8; ++i)
    obs << double(i);
  return obs;
}

} // namespace

TEST(AleaRegressions, RebinIntoOneBinKeepsTheSum)
{
  alps::RealObservable obs = one_to_eight();
  ASSERT_EQ(obs.bin_number(), 8u);
  obs.set_bin_number(1);
  ASSERT_EQ(obs.bin_number(), 1u);
  EXPECT_DOUBLE_EQ(obs.bin_value(0), 36.);
}

TEST(AleaRegressions, RebinWithAndWithoutPartialBin)
{
  alps::RealObservable full = one_to_eight();
  full.set_bin_number(2);
  ASSERT_EQ(full.bin_number(), 2u);
  EXPECT_DOUBLE_EQ(full.bin_value(0), 10.);
  EXPECT_DOUBLE_EQ(full.bin_value(1), 26.);

  alps::RealObservable partial = one_to_eight();
  partial.set_bin_number(3);
  // bins of 3: {1,2,3}, {4,5,6} and the partly filled {7,8}, which
  // bin_number() does not count
  ASSERT_EQ(partial.bin_number(), 2u);
  EXPECT_EQ(partial.bin_size(), 3u);
  EXPECT_DOUBLE_EQ(partial.bin_value(0), 6.);
  EXPECT_DOUBLE_EQ(partial.bin_value(1), 15.);
}

TEST(AleaRegressions, McdataAssignmentCopiesMaxBinNumber)
{
  alps::RealObservable obs = one_to_eight();
  obs.set_bin_number(4);
  alps::alea::mcdata<double> a(obs);
  ASSERT_EQ(a.max_bin_number(), 4u);

  alps::alea::mcdata<double> b;
  ASSERT_EQ(b.max_bin_number(), 0u);
  b = a;
  EXPECT_EQ(b.max_bin_number(), a.max_bin_number());
}

TEST(AleaRegressions, IntegerHistogramCoversWholeRange)
{
  alps::IntHistogramObservable h("h", 0, 5, 2);
  ASSERT_EQ(h.size(), 3u);
  h << 4;
  EXPECT_EQ(h[2], 1u);
  EXPECT_EQ(h.count(), 1u);
}

TEST(AleaRegressions, FloatingHistogramDropsValuesPastLastBin)
{
  alps::RealHistogramObservable h("h", 0., 2.3, 1.);
  ASSERT_EQ(h.size(), 2u);
  h << 2.2;
  EXPECT_EQ(h.count(), 0u);
}

TEST(AleaRegressions, HistogramXmlRejectsEntriesBeyondNvalues)
{
  std::istringstream valid("<HISTOGRAM nvalues=\"2\">"
                           "<ENTRY><COUNT>3</COUNT><VALUE>1</VALUE></ENTRY>"
                           "<ENTRY><COUNT>3</COUNT><VALUE>2</VALUE></ENTRY>"
                           "</HISTOGRAM>");
  alps::XMLTag tag = alps::parse_tag(valid);
  alps::HistogramObservableData<alps::int32_t> data(valid, tag);
  ASSERT_EQ(data.size(), 2u);
  EXPECT_EQ(data[1], 2u);

  std::istringstream overlong("<HISTOGRAM nvalues=\"1\">"
                              "<ENTRY><COUNT>3</COUNT><VALUE>1</VALUE></ENTRY>"
                              "<ENTRY><COUNT>3</COUNT><VALUE>2</VALUE></ENTRY>"
                              "</HISTOGRAM>");
  tag = alps::parse_tag(overlong);
  EXPECT_THROW((alps::HistogramObservableData<alps::int32_t>(overlong, tag)), std::runtime_error);
}
