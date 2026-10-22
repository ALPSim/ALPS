#include <gtest/gtest.h>
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include "dmtk/dmtk.h"
#include <cstring>
#include <new>

namespace {

typedef dmtk::System<double> System;

struct Calls {
  int count;
  bool result;
};

bool record(System&, size_t, void* data)
{
  Calls& calls = *static_cast<Calls*>(data);
  ++calls.count;
  return calls.result;
}

} // namespace

TEST(DmrgSignals, DefaultConstructedSignalIsNull)
{
  // Build in poisoned storage so members the constructor skips stay nonzero.
  alignas(dmtk::Signal<System>) unsigned char storage[sizeof(dmtk::Signal<System>)];
  std::memset(storage, 0xAB, sizeof(storage));
  dmtk::Signal<System>* signal = new (storage) dmtk::Signal<System>();
  EXPECT_EQ(signal->data(), nullptr);
  EXPECT_EQ(signal->signal_id(), 0u);
  signal->~Signal();
}

TEST(DmrgSignals, UnknownNameConnectsNothing)
{
  System S;
  Calls calls = {0, false};
  S.signal_connect(record, std::string("NO_SUCH_SIGNAL"), &calls);
  // get_signal_id() returns -1 here; it must not become a size_t id.
  EXPECT_TRUE(S.signal_emit(static_cast<size_t>(-1)));
  EXPECT_EQ(calls.count, 0);
}

TEST(DmrgSignals, EmitByName)
{
  System S;
  Calls calls = {0, false};
  S.signal_connect(record, std::string("SYSTEM_SIGNAL_GS"), &calls);
  EXPECT_EQ(S.signal_id("SYSTEM_SIGNAL_GS"), dmtk::SYSTEM_SIGNAL_GS);
  EXPECT_FALSE(S.signal_emit(std::string("SYSTEM_SIGNAL_GS")));
  EXPECT_EQ(calls.count, 1);
  // An unknown name runs no handler and vetoes nothing.
  EXPECT_TRUE(S.signal_emit(std::string("NO_SUCH_SIGNAL")));
  EXPECT_EQ(calls.count, 1);
}
