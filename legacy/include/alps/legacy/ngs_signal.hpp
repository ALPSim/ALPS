// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Deprecated alps::ngs name, included by <alps/utilities/signal.hpp> when the
// legacy layer is installed.
#pragma once

namespace alps {
    namespace ngs {
        using signal [[deprecated("use alps::signal from <alps/utilities/signal.hpp>")]] = ::alps::signal;
    }
}
