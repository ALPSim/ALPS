/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/temperature_scan.h>
#include <gtest/gtest.h>
#include <vector>

class RecordingWorker {
public:
  inline static std::vector<double> temperatures;
  RecordingWorker(const alps::Parameters &) {}
  void init_observables(const alps::Parameters &, const alps::ObservableSet &) {}
  void run(const alps::ObservableSet &) {}
  void set_beta(double beta) { temperatures.push_back(1 / beta); }
  void save(alps::ODump &) const {}
  void load(alps::IDump &) {}
};
TEST(TemperatureScan, UsesInitialThermalizationThenVisitsEveryTemperature) {
  alps::Parameters params;
  params["THERMALIZATION"] = 3;
  params["SWEEPS"] = 5;
  params["NUM_TEMPERATURES"] = 5;
  params["INITIAL_TEMPERATURE"] = .1;
  params["DIFF_TEMPERATURE"] = .1;
  params["INITIAL_THERMALIZATION"] = 10;
  RecordingWorker::temperatures.clear();
  std::vector<alps::ObservableSet> observables;
  alps::parapack::temperature_scan_adaptor<RecordingWorker> worker(params);
  worker.init_observables(params, observables);
  ASSERT_EQ(observables.size(), 5);
  EXPECT_DOUBLE_EQ(worker.progress(), 0);
  for (int step = 0; step < 47; ++step) {
    ASSERT_LT(worker.progress(), 1);
    worker.run(observables);
    const double expected = step < 15 ? (step + 1) / 75. : .2 + (step - 14) / 40.;
    EXPECT_NEAR(worker.progress(), expected, 1e-14) << "step " << step;
  }
  EXPECT_DOUBLE_EQ(worker.progress(), 1);
  ASSERT_EQ(RecordingWorker::temperatures.size(), 47);
  for (int step = 0; step < 47; ++step) {
    const double expected = step < 15 ? .1 : .2 + ((step - 15) / 8) * .1;
    EXPECT_NEAR(RecordingWorker::temperatures[step], expected, 1e-14) << "step " << step;
  }
  worker.run(observables);
  EXPECT_EQ(RecordingWorker::temperatures.size(), 47);
}
