// SPDX-License-Identifier: MIT
//
// Parapack's MPI messages: the in-memory osiris dump in <alps/parapack/mpi.h>
// must give back every field it was given, and must reject a message whose
// length does not match what load() reads.

#include <alps/parapack/clone.h>
#include <alps/parapack/mpi.h>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

alps::Parameters example_parameters() {
  alps::Parameters params;
  params["LATTICE"] = "square lattice";
  params["L"] = 16;
  params["T"] = 0.125;
  params["SEED"] = 4711;
  params["DISORDER_SEED"] = 815;
  params["MODEL"] = "spin";
  return params;
}

alps::clone_info example_clone_info() {
  // Initializing sets the checkpoint name and the seeds from SEED.
  alps::clone_info info(3, example_parameters(), "run/task7");
  info.start("equilibrating");
  info.stop();
  info.start("running");
  info.stop();
  info.set_progress(0.375);
  return info;
}

template<typename T>
std::string text(T const& value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

void expect_equal(alps::clone_create_msg_t const& a, alps::clone_create_msg_t const& b) {
  EXPECT_EQ(a.task_id, b.task_id);
  EXPECT_EQ(a.clone_id, b.clone_id);
  EXPECT_EQ(a.group_id, b.group_id);
  EXPECT_EQ(text(a.params), text(b.params));
  EXPECT_EQ(a.base, b.base);
  EXPECT_EQ(a.is_new, b.is_new);
}

void expect_equal(alps::clone_info_msg_t const& a, alps::clone_info_msg_t const& b) {
  EXPECT_EQ(a.task_id, b.task_id);
  EXPECT_EQ(a.clone_id, b.clone_id);
  EXPECT_EQ(a.group_id, b.group_id);
  EXPECT_EQ(a.info.clone_id(), b.info.clone_id());
  EXPECT_EQ(a.info.progress(), b.info.progress());
  EXPECT_EQ(a.info.checkpoints(), b.info.checkpoints());
  EXPECT_EQ(a.info.worker_seed(), b.info.worker_seed());
  EXPECT_EQ(a.info.disorder_seed(), b.info.disorder_seed());
  // The clone's own host list is not part of the message, as it was not
  // with Boost serialization, so compare the phases rather than write_xml().
  ASSERT_EQ(a.info.phases().size(), b.info.phases().size());
  for (std::size_t i = 0; i < a.info.phases().size(); ++i) {
    alps::clone_phase const& p = a.info.phases()[i];
    alps::clone_phase const& q = b.info.phases()[i];
    EXPECT_EQ(p.hosts(), q.hosts());
    EXPECT_EQ(p.user(), q.user());
    EXPECT_EQ(p.phase(), q.phase());
    EXPECT_EQ(p.start_time(), q.start_time());
    EXPECT_EQ(p.stop_time(), q.stop_time());
  }
}

template<typename T>
T round_trip(T const& value) {
  T result;
  alps::mpi::detail::deserialize(alps::mpi::detail::serialize(&value, 1), &result, 1);
  return result;
}

// save() writes two values, load() reads only the first
struct reads_less {
  int a = 0, b = 0;
  void save(alps::ODump& dp) const { dp << a << b; }
  void load(alps::IDump& dp) { dp >> a; }
};

// save() writes one value, load() reads two
struct reads_more {
  int a = 0, b = 0;
  void save(alps::ODump& dp) const { dp << a; }
  void load(alps::IDump& dp) { dp >> a >> b; }
};

} // namespace

TEST(ParapackMessages, CloneCreateRoundTrip) {
  alps::clone_create_msg_t const msg(7, 3, 2, example_parameters(), "run/task7", true);
  expect_equal(round_trip(msg), msg);
}

TEST(ParapackMessages, CloneInfoRoundTrip) {
  alps::clone_info_msg_t const msg(7, 3, 2, example_clone_info());
  ASSERT_EQ(msg.info.phases().size(), 2u);
  expect_equal(round_trip(msg), msg);
}

TEST(ParapackMessages, CloneHaltRoundTrip) {
  alps::clone_halt_msg_t const msg(7, 3, 2);
  alps::clone_halt_msg_t const result = round_trip(msg);
  EXPECT_EQ(result.task_id, msg.task_id);
  EXPECT_EQ(result.clone_id, msg.clone_id);
  EXPECT_EQ(result.group_id, msg.group_id);
}

TEST(ParapackMessages, ParametersRoundTrip) {
  alps::Parameters const params = example_parameters();
  EXPECT_EQ(text(round_trip(params)), text(params));
}

TEST(ParapackMessages, StringsRoundTrip) {
  // One string goes as its characters, several through the dump.
  std::string const one("a single string");
  EXPECT_EQ(round_trip(one), one);
  std::vector<std::string> const many{"first", "", "third string"};
  std::vector<std::string> result(many.size());
  alps::mpi::detail::deserialize(alps::mpi::detail::serialize(many.data(), 3),
                                 result.data(), 3);
  EXPECT_EQ(result, many);
}

TEST(ParapackMessages, BoolArrayRoundTrip) {
  // bool is not sent as an MPI type, so it goes through the dump.
  bool const values[] = {true, false, true, true};
  bool result[4] = {};
  alps::mpi::detail::deserialize(alps::mpi::detail::serialize(values, 4), result, 4);
  EXPECT_TRUE(std::equal(values, values + 4, result));
}

TEST(ParapackMessages, RejectsTruncatedMessage) {
  alps::clone_create_msg_t const msg(7, 3, 2, example_parameters(), "run/task7", true);
  std::string buffer = alps::mpi::detail::serialize(&msg, 1);
  buffer.pop_back();
  alps::clone_create_msg_t result;
  EXPECT_THROW(alps::mpi::detail::deserialize(buffer, &result, 1), std::runtime_error);
}

TEST(ParapackMessages, RejectsLoadThatReadsLessThanSaveWrote) {
  reads_less const value{1, 2};
  reads_less result;
  std::string const buffer = alps::mpi::detail::serialize(&value, 1);
  EXPECT_THROW(alps::mpi::detail::deserialize(buffer, &result, 1), std::runtime_error);
}

TEST(ParapackMessages, RejectsLoadThatReadsMoreThanSaveWrote) {
  reads_more const value{1, 2};
  reads_more result;
  std::string const buffer = alps::mpi::detail::serialize(&value, 1);
  EXPECT_THROW(alps::mpi::detail::deserialize(buffer, &result, 1), std::runtime_error);
}

TEST(ParapackMessages, SendsMessagesBetweenRanks) {
  alps::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  int const tag = 17;
  alps::clone_create_msg_t const create(7, 3, 2, example_parameters(), "run/task7", true);
  alps::clone_info_msg_t const info(7, 3, 2, example_clone_info());
  if (world.rank() == 0) {
    for (int p = 1; p < world.size(); ++p) {
      alps::mpi::send(world, p, tag, create);
      alps::mpi::send(world, p, tag + 1, info);
    }
  } else {
    alps::clone_create_msg_t create_received;
    alps::mpi::status const st = alps::mpi::recv(world, 0, tag, create_received);
    EXPECT_EQ(st.source(), 0);
    EXPECT_EQ(st.tag(), tag);
    expect_equal(create_received, create);
    alps::clone_info_msg_t info_received;
    alps::mpi::recv(world, 0, tag + 1, info_received);
    expect_equal(info_received, info);
  }
  world.barrier();
}

TEST(ParapackMessages, BroadcastsParameters) {
  alps::mpi::communicator world;
  alps::Parameters params;
  if (world.rank() == 0) params = example_parameters();
  alps::mpi::broadcast_object(world, params, 0);
  EXPECT_EQ(text(params), text(example_parameters()));
}
