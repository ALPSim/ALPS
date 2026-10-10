// SPDX-License-Identifier: MIT
#include <alps/utility/citations.hpp>
#include <alps/utility/copyright.hpp>
#include <cstdlib>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>

TEST(Citations, NoticesAndSuppression) {
  auto set_no_citations = [](const char* value) {
#ifdef _WIN32
    _putenv_s("ALPS_NO_CITATIONS", value ? value : "");
#else
    if (value) setenv("ALPS_NO_CITATIONS", value, 1);
    else unsetenv("ALPS_NO_CITATIONS");
#endif
  };
  const std::string title = "The ALPS project release 3.0:";
  for (const char* component : {"framework", "scheduler", "dmrg", "qwl", "fulldiag",
                               "sparsediag", "spinmc", "worm", "dirloop_sse", "looper",
                               "dmft", "interaction", "hybridization", "hirschfye"}) {
    const std::string notice = alps::citation_text(component);
    const std::string details = alps::citation_details(component);
    const auto first = notice.find(title);
    if (first == std::string::npos || notice.find(title, first + title.size()) != std::string::npos)
      throw std::runtime_error(std::string(component) + ": expected one framework reference");
    if (notice.find("P05001") != std::string::npos)
      throw std::runtime_error("Outdated framework paper in notice");
    set_no_citations(nullptr);
    std::ostringstream enabled;
    alps::print_copyright(enabled, component);
    const std::string copyright = enabled.str().substr(0, enabled.str().size() - notice.size());
    if (copyright.find("copyright (c)") == std::string::npos ||
        copyright.find("Licensed under the MIT License.") == std::string::npos ||
        copyright.find("https://github.com/ALPSim/ALPS/blob/master/LICENSE.txt") == std::string::npos)
      throw std::runtime_error("Copyright/license banner missing");
    for (const char* value : {static_cast<const char*>(nullptr), "", "0", "1", "true", "01", "1 "}) {
      set_no_citations(value);
      const std::string expected = value && std::string(value) == "1" ? "" : notice;
      std::ostringstream compact, full, banner;
      alps::print_citations(compact, component);
      alps::print_citation_details(full, component);
      alps::print_copyright(banner, component);
      if (compact.str() != expected || full.str() != details ||
          alps::citation_text(component) != notice || alps::citation_details(component) != details)
        throw std::runtime_error("ALPS_NO_CITATIONS changed the wrong citation output");
      if (banner.str() != copyright + expected)
        throw std::runtime_error("ALPS_NO_CITATIONS did not preserve the copyright/license banner");
    }
  }
  set_no_citations(nullptr);
  if (alps::citation_details("interaction").find("10.1103/PhysRevB.72.035122") == std::string::npos)
    throw std::runtime_error("CT-INT algorithm reference missing");
  if (alps::citation_details("hybridization").find("10.1103/PhysRevB.72.035122") != std::string::npos)
    throw std::runtime_error("CT-HYB incorrectly recommends CT-INT");
  EXPECT_THROW(alps::citation_text("unknown-component"), std::invalid_argument);
}
