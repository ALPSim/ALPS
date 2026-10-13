#include "observable_checks.hpp"
#include <string>
#include <vector>
#include <alps/alea.h>
#include <boost/random.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

TEST(AleaMerge, divide) {
  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int(1u); // Preserve the historical default seed.
  random_type random(random_int);

  const int MCS = 128;
  const int nsets = 2;

  std::vector<alps::ObservableSet> obssets(nsets);
  BOOST_FOREACH(alps::ObservableSet &obs, obssets){
    obs << alps::RealObservable("two");
    obs.reset(true);
    for(int i = 0; i < MCS; ++i) obs["two"] << 2.0 + (2 * random() - 1) * 0.5;
  }

  // merge
  for(int i = 1; i < nsets; ++i) obssets[0] << obssets[i];
  alps::ObservableSet& obs = obssets[0];

  alps::RealObsevaluator two = obs["two"];
  alps::RealObsevaluator result("4.0 / two");
  result = 4.0 / two;
  obs.addObservable(result);
  alps_test::expect_estimate(result, {{1.98611, 0.0182, 5.00001e-06, 5.00001e-05}}, "4.0 / two");
  EXPECT_EQ(result.count(), 256u);
  EXPECT_EQ(result.converged_errors(), alps::MAYBE_CONVERGED);
  alps_test::expect_estimate(two, {{2.01382, 0.0184, 5.00001e-06, 5.00001e-05}}, "two");
  EXPECT_EQ(two.count(), 256u);
  EXPECT_EQ(two.converged_errors(), alps::MAYBE_CONVERGED);
  EXPECT_NEAR(two.tau(), 0.022, 0.000500001);


}
