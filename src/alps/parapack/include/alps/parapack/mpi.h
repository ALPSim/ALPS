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

// Point-to-point messages and the few collectives that parapack needs, on
// top of alps::mpi::communicator from <alps/utilities/mpi.hpp>. They replace
// the Boost.MPI calls parapack used before:
//
//   comm.send(dest, tag)          ->  mpi::send(comm, dest, tag)
//   comm.send(dest, tag, value)   ->  mpi::send(comm, dest, tag, value)
//   comm.recv(src, tag, value)    ->  mpi::recv(comm, src, tag, value)
//   comm.iprobe(src, tag)         ->  mpi::iprobe(comm, src, tag)
//   comm.split(color)             ->  mpi::split(comm, color)
//   broadcast(comm, object, root) ->  mpi::broadcast_object(comm, object, root)
//   gather(comm, string, ...)     ->  mpi::gather(comm, string, ...)
//
// Arithmetic values and arrays go as their MPI type. Other values are
// written to a Boost.Serialization binary archive and sent as one MPI_BYTE
// message, so their type needs a serialize() function, as with Boost.MPI.

#ifndef PARAPACK_MPI_H
#define PARAPACK_MPI_H

#include <alps/config.h>

#ifdef ALPS_HAVE_MPI

#include <alps/utilities/mpi.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace alps {
namespace mpi {

int const any_source = MPI_ANY_SOURCE;
int const any_tag = MPI_ANY_TAG;

/// Source and tag of a received or probed message
class status {
public:
  explicit status(MPI_Status const& s) : status_(s) {}
  int source() const { return status_.MPI_SOURCE; }
  int tag() const { return status_.MPI_TAG; }
  operator MPI_Status const&() const { return status_; }
private:
  MPI_Status status_;
};

/// Returns the status of a matching message if one has arrived
inline std::optional<status> iprobe(communicator const& comm, int source, int tag) {
  int flag;
  MPI_Status s;
  MPI_Iprobe(source, tag, comm, &flag, &s);
  if (flag) return status(s);
  return std::nullopt;
}

/// Waits for a matching message and returns its status
inline status probe(communicator const& comm, int source, int tag) {
  MPI_Status s;
  MPI_Probe(source, tag, comm, &s);
  return status(s);
}

/// Splits the communicator by color, keeping the rank order. Processes that
/// pass MPI_UNDEFINED get a communicator that compares equal to MPI_COMM_NULL.
inline communicator split(communicator const& comm, int color, int key) {
  MPI_Comm c;
  MPI_Comm_split(comm, color, key, &c);
  if (c == MPI_COMM_NULL) return communicator(c, comm_attach);
  return communicator(c, take_ownership);
}

inline communicator split(communicator const& comm, int color) {
  return split(comm, color, comm.rank());
}

namespace detail {

// Values sent with their own MPI type; bool has none in every MPI.
template<typename T>
inline constexpr bool is_mpi_value = std::is_arithmetic<T>::value && !std::is_same<T, bool>::value;

// A byte count or offset as the int that MPI takes. The message goes to stderr
// too, since parapack's MPI programs abort while unwinding, before their catch
// blocks can print it (see alps::mpi::environment).
inline int mpi_count(std::size_t n) {
  if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    char const* const message = "alps::mpi: message exceeds the MPI count limit of 2 GiB";
    std::cerr << "Error: " << message << std::endl;
    throw std::length_error(message);
  }
  return static_cast<int>(n);
}

template<typename T>
std::string serialize(T const* values, int n) {
  std::ostringstream os;
  {
    boost::archive::binary_oarchive ar(os);
    for (int i = 0; i < n; ++i) ar << values[i];
  }
  return os.str();
}

template<typename T>
void deserialize(std::string const& buffer, T* values, int n) {
  std::istringstream is(buffer);
  boost::archive::binary_iarchive ar(is);
  for (int i = 0; i < n; ++i) ar >> values[i];
}

} // namespace detail

/// Sends a message without data
inline void send(communicator const& comm, int dest, int tag) {
  MPI_Send(nullptr, 0, MPI_BYTE, dest, tag, comm);
}

/// Receives a message without data
inline status recv(communicator const& comm, int source, int tag) {
  MPI_Status s;
  MPI_Recv(nullptr, 0, MPI_BYTE, source, tag, comm, &s);
  return status(s);
}

/// Sends n values: arithmetic values as themselves, anything else serialized
template<typename T>
void send(communicator const& comm, int dest, int tag, T const* values, int n) {
  if constexpr (detail::is_mpi_value<T>) {
    MPI_Send(const_cast<T*>(values), n, detail::mpi_type<T>(), dest, tag, comm);
  } else {
    std::string const buffer = detail::serialize(values, n);
    MPI_Send(const_cast<char*>(buffer.data()), detail::mpi_count(buffer.size()), MPI_BYTE,
             dest, tag, comm);
  }
}

/// Receives n values sent with send(comm, dest, tag, values, n)
template<typename T>
status recv(communicator const& comm, int source, int tag, T* values, int n) {
  if constexpr (detail::is_mpi_value<T>) {
    MPI_Status s;
    MPI_Recv(values, n, detail::mpi_type<T>(), source, tag, comm, &s);
    return status(s);
  } else {
    // Probe first, so that the receive matches the message whose size we read.
    status const st = probe(comm, source, tag);
    int count;
    MPI_Get_count(&static_cast<MPI_Status const&>(st), MPI_BYTE, &count);
    std::string buffer(count, '\0');
    MPI_Recv(&buffer[0], count, MPI_BYTE, st.source(), st.tag(), comm, MPI_STATUS_IGNORE);
    detail::deserialize(buffer, values, n);
    return st;
  }
}

/// Sends one value
template<typename T>
void send(communicator const& comm, int dest, int tag, T const& value) {
  send(comm, dest, tag, &value, 1);
}

/// Receives one value sent with send(comm, dest, tag, value)
template<typename T>
status recv(communicator const& comm, int source, int tag, T& value) {
  return recv(comm, source, tag, &value, 1);
}

/// Broadcasts a serializable object from root
template<typename T>
void broadcast_object(communicator const& comm, T& value, int root) {
  std::string buffer;
  if (comm.rank() == root) buffer = detail::serialize(&value, 1);
  broadcast(comm, buffer, root);
  if (comm.rank() != root) detail::deserialize(buffer, &value, 1);
}

/// Gathers one string from each process into values at root; values is
/// left alone on the other processes
inline void gather(communicator const& comm, std::string const& value,
                   std::vector<std::string>& values, int root) {
  int const n = comm.size();
  int size = detail::mpi_count(value.size());
  std::vector<int> sizes(n), offsets(n);
  MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, root, comm);
  std::size_t total = 0;
  for (int p = 0; p < n; ++p) {
    offsets[p] = detail::mpi_count(total);
    total += sizes[p];
  }
  std::vector<char> buffer(total);
  MPI_Gatherv(const_cast<char*>(value.data()), size, MPI_CHAR, buffer.data(), sizes.data(),
              offsets.data(), MPI_CHAR, root, comm);
  if (comm.rank() != root) return;
  values.resize(n);
  for (int p = 0; p < n; ++p) values[p].assign(buffer.data() + offsets[p], sizes[p]);
}

/// Sends one string to root, which calls the other gather()
inline void gather(communicator const& comm, std::string const& value, int root) {
  int size = detail::mpi_count(value.size());
  MPI_Gather(&size, 1, MPI_INT, nullptr, 1, MPI_INT, root, comm);
  MPI_Gatherv(const_cast<char*>(value.data()), size, MPI_CHAR, nullptr, nullptr, nullptr,
              MPI_CHAR, root, comm);
}

} // namespace mpi
} // namespace alps

#endif // ALPS_HAVE_MPI

#endif // PARAPACK_MPI_H
