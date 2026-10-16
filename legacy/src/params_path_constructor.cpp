// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// The ALPS 2.x text-file constructor, kept as a thin layer over the explicit
// parser so the two cannot diverge.
#include <alps/ngs/params_from_file.hpp>

namespace alps {

    params::params(boost::filesystem::path const & path)
        : params(params_from_file(path))
    {}

}
