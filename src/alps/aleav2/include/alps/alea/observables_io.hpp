/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

/**
 * Observables in files and messages. See SCHEMA.md for the layouts.
 *
 * There are two kinds of files. Results are for programs that evaluate a
 * simulation: mean, error, autocorrelation time and bins, in a layout that
 * stays readable. Checkpoints are for the simulation itself, to continue
 * where it stopped, and are read only by the version that wrote them.
 */
#pragma once

#include <alps/alea/observables.hpp>

#include <string>

namespace alps { namespace hdf5 { class archive; } }

namespace alps { namespace alea {

/** Write results below `path`, one group per observable */
void save_results(hdf5::archive &ar, const std::string &path, const result_set &results);

/**
 * Read results below `path`: those written by `save_results`, and the result
 * files of ALPS 2 and ALPS 3.0. Means and errors are returned as stored.
 * Results with stored bins can enter functions of results. Results read from
 * a file cannot be merged with other runs.
 */
result_set load_results(hdf5::archive &ar, const std::string &path);

/**
 * One result in the group `group`, as `save_results` writes it. Of a signed
 * observable this keeps mean and error; its samples are written with the set.
 */
void save_result(hdf5::archive &ar, const std::string &group, const observable_result &result);
observable_result load_result(hdf5::archive &ar, const std::string &group);

/** Write everything needed to continue measuring */
void save_checkpoint(hdf5::archive &ar, const std::string &path, const observable_set &set);

/**
 * Continue the observables saved with `save_checkpoint`. They replace those of
 * the same name in `set`. Anything else below `path`, such as the measurements
 * of an earlier version of ALPS, is refused.
 */
void load_checkpoint(hdf5::archive &ar, const std::string &path, observable_set &set);

/**
 * Results as bytes, to send them to another process of the same program,
 * where `unpack` returns them. Unlike results read from a file they can be
 * merged there.
 */
std::string pack(const result_set &results);
result_set unpack(const std::string &bytes);

}} /* namespace alps::alea */
