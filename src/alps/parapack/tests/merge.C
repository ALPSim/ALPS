#include <alps/parapack/measurement.h>
#include <gtest/gtest.h>
#include <cmath>

// The numbers 0 to 99, each measured 100 times by its own clone. As a
// disorder average every clone is one sample: mean 49.5 and error
// sqrt(833.25 / 99) = 2.9011, from the 100 uncorrelated clone means.
// Merged with equal weight per sample all 10000 samples count.
TEST(MeasurementMerge, DisorderAverageTakesOneSamplePerClone) {
  alps::alea::observable_set merged, random_clone;
  for (int i = 0; i < 100; ++i) {
    alps::alea::observable_set clone;
    clone << alps::alea::observable("obs");
    for (int j = 0; j < 100; ++j)
      clone["obs"] << static_cast<double>(i);
    alps::merge_clone(merged, clone.results(), false);
    alps::merge_random_clone(random_clone, clone.results());
  }
  const alps::alea::result_set disorder = random_clone.results();
  EXPECT_EQ(disorder["obs"].count(), 100u);
  EXPECT_DOUBLE_EQ(disorder["obs"].mean(), 49.5);
  EXPECT_NEAR(disorder["obs"].error(), std::sqrt(833.25 / 99), 1e-12);

  const alps::alea::result_set all = merged.results();
  EXPECT_EQ(all["obs"].count(), 10000u);
  EXPECT_NEAR(all["obs"].mean(), 49.5, 1e-12);
}
