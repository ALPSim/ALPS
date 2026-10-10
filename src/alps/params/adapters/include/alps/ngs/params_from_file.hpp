// SPDX-License-Identifier: MIT
#ifndef ALPS_NGS_PARAMS_FROM_FILE_HPP
#define ALPS_NGS_PARAMS_FROM_FILE_HPP

#include <alps/ngs/params.hpp>
#include <boost/filesystem/path.hpp>

namespace alps {
// Parse the legacy text grammar. Requires ALPS::alps.
ALPS_DECL params params_from_file(boost::filesystem::path const& filename);
}
#endif
