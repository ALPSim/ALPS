// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Deprecated macro names, included by <alps/utilities/stringify.hpp> when the
// legacy layer is installed. Macros cannot carry [[deprecated]]; these
// disappear with the legacy layer.
#pragma once

#define ALPS_NGS_STRINGIFY(arg) ALPS_STRINGIFY(arg)
#define ALPS_NGS_STRINGIFY_HELPER(arg) ALPS_STRINGIFY_HELPER(arg)
