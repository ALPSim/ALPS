// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Deprecated alps::params member, included inside the class body by
// <alps/ngs/params.hpp> when the legacy layer is installed. Defined in
// legacy/src/params_path_constructor.cpp, part of ALPS::alps.

[[deprecated("use alps::params_from_file(path) from <alps/ngs/params_from_file.hpp>; link ALPS::alps")]]
ALPS_DECL params(boost::filesystem::path const & path);
