#ifndef PARAPACK_STAGING_H
#define PARAPACK_STAGING_H

#include <filesystem>
#include <boost/tuple/tuple.hpp>
#include <queue>
#include <alps/export.h>

namespace alps {
namespace parapack {

typedef
boost::tuple<int /* task_id */,int /* clone_id */,int /* group_id*/> suspended_queue_t;

ALPS_DECL void load_checkpoints(std::filesystem::path const& file_chp,
		      std::filesystem::path const& basedir,
		      std::queue<suspended_queue_t>& suspended_queue);

} // end namespace parapack
} // end namespace alps
#endif // PARAPACK_STAGING_H
