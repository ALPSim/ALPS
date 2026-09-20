#include <boost/filesystem/fstream.hpp>
#include <stdexcept>
#include "staging.h"

#include <string>
#include <fstream>
#include <iostream>

namespace alps {
namespace parapack {

void load_checkpoints(boost::filesystem::path const& file_chp,
		      boost::filesystem::path const& basedir,
		      std::queue<suspended_queue_t>& suspended_queue) {

  std::cout << "  tasks ordered by " << file_chp << " = ";
  boost::filesystem::ifstream stream(file_chp);
  if (!stream) {
    std::cout << "no" << std::endl;
  } else {
    std::cout << "yes" << std::endl;
    int task, clone, group;
    while (stream >> task >> clone >> group)
      suspended_queue.push(boost::tuple<int,int,int>(task, clone, group));
    if (!stream.eof())
      throw std::runtime_error("Invalid checkpoint staging entry: " + file_chp.string());
  }
}

}
}

