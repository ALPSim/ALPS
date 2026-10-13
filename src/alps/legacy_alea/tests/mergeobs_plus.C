#include "observable_checks.hpp"
#include <string>
#include <vector>
#include <alps/alea.h>
#include <boost/random.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

TEST(AleaMerge, plus) {
  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int(1u); // Preserve the historical default seed.
  random_type random(random_int);

  const int MCS = 128;
  const int nsets = 2;

  std::vector<alps::ObservableSet> obssets(nsets);
  BOOST_FOREACH(alps::ObservableSet &obs, obssets){
    obs << alps::RealObservable("one");
    obs.reset(true);
    for(int i = 0; i < MCS; ++i) obs["one"] << 1.0 + (2 * random() - 1) * 0.3;
  }

  // merge
  for(int i = 1; i < nsets; ++i) obssets[0] << obssets[i];
  alps::ObservableSet& obs = obssets[0];

  alps::RealObsevaluator one = obs["one"];
  alps::RealObsevaluator result("one + 1.0");
  result = one + 1.0;
  obs.addObservable(result);
  alps_test::expect_estimate(one, {{1.00829, 0.011, 5.00001e-06, 0.000500001}}, "one");
  EXPECT_EQ(one.count(), 256u);
  EXPECT_EQ(one.converged_errors(), alps::MAYBE_CONVERGED);
  EXPECT_NEAR(one.tau(), 0.022, 0.000500001);
  alps_test::expect_estimate(result, {{2.00829, 0.011, 5.00001e-06, 0.000500001}}, "one + 1.0");
  EXPECT_EQ(result.count(), 256u);
  EXPECT_EQ(result.converged_errors(), alps::MAYBE_CONVERGED);
  EXPECT_NEAR(result.tau(), -0.471, 0.000500001);


}
