// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Deprecated alps::ngs name, included by <alps/utilities/stacktrace.hpp> when
// the legacy layer is installed.
#pragma once

#include <string>

namespace alps {
    namespace ngs {
        [[deprecated("use alps::stacktrace() from <alps/utilities/stacktrace.hpp>")]]
        inline std::string stacktrace() { return ::alps::stacktrace(); }
    }
}
