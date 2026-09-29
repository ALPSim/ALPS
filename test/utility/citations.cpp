// SPDX-License-Identifier: MIT
#include <alps/utility/citations.hpp>
#include <alps/utility/copyright.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

int main() {
  const std::string title = "The ALPS project release 3.0:";
  for (const char* component : {"framework", "scheduler", "dmrg", "qwl", "fulldiag",
                               "sparsediag", "spinmc", "worm", "dirloop_sse", "looper",
                               "dmft", "interaction", "hybridization", "hirschfye"}) {
    std::ostringstream out;
    alps::print_copyright(out, component);
    const std::string text = out.str();
    const auto first = text.find(title);
    if (first == std::string::npos || text.find(title, first + title.size()) != std::string::npos)
      throw std::runtime_error(std::string(component) + ": expected one framework reference");
    if (text.find("P05001") != std::string::npos)
      throw std::runtime_error("Outdated framework paper in notice");
  }
  if (alps::citation_details("interaction").find("10.1103/PhysRevB.72.035122") == std::string::npos)
    throw std::runtime_error("CT-INT algorithm reference missing");
  if (alps::citation_details("hybridization").find("10.1103/PhysRevB.72.035122") != std::string::npos)
    throw std::runtime_error("CT-HYB incorrectly recommends CT-INT");
  try {
    alps::citation_text("unknown-component");
  } catch (const std::invalid_argument&) {
    return 0;
  }
  throw std::runtime_error("Unknown citation component was silently accepted");
}
