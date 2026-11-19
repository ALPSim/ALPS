/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

/**
 * Header written with the state of a running accumulator.
 *
 * The state is what is needed to continue an accumulator, bit for bit, where
 * it was saved. It is not a result, and its layout may change between
 * versions: a state is read only by code with the same `STATE_VERSION`.
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/serialize.hpp>

#include <stdexcept>
#include <string>

namespace alps { namespace alea { namespace internal {

const uint64_t STATE_VERSION = 1;

enum state_type {
    STATE_MEAN = 1, STATE_VAR = 2, STATE_COV = 3, STATE_AUTOCORR = 4, STATE_BATCH = 5
};

inline void write_state_header(serializer &s, state_type type)
{
    serialize(s, "@alea_state", STATE_VERSION);
    serialize(s, "@alea_type", static_cast<uint64_t>(type));
}

inline void read_state_header(deserializer &s, state_type type)
{
    uint64_t version, stored_type;
    deserialize(s, "@alea_state", version);
    if (version != STATE_VERSION) {
        throw std::runtime_error(
            "Accumulator state has version " + std::to_string(version)
            + ", this code reads version " + std::to_string(STATE_VERSION)
            + ". A state can only be continued by the version that wrote it.");
    }
    deserialize(s, "@alea_type", stored_type);
    if (stored_type != static_cast<uint64_t>(type)) {
        throw std::runtime_error(
            "Accumulator state is of type " + std::to_string(stored_type)
            + ", expected type " + std::to_string(static_cast<uint64_t>(type)));
    }
}

inline uint64_t read_uint(deserializer &s, const std::string &key)
{
    uint64_t value;
    deserialize(s, key, value);
    return value;
}

}}} /* namespace alps::alea::internal */
