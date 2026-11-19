/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_SIMULATION_P_H
#define PARAPACK_SIMULATION_P_H

#include "clone_info_p.h"

#include <alps/parameter/parameters_p.h>
#include <alps/alea/observables.hpp>
#include <alps/scheduler/montecarlo.h>
#include <alps/parser/xmlhandler.h>
#include <filesystem>

namespace alps {

struct simulation_xml_writer {
  simulation_xml_writer(std::filesystem::path file, bool write_xml, bool make_backup,
    Parameters const& params, std::vector<alps::alea::result_set> const& obs,
    std::deque<clone_info> const& info) {
    std::filesystem::path file_bak(file.parent_path() / (file.filename().string() + ".bak"));
    if (make_backup && exists(file)) rename(file, file_bak);
    oxstream os(file);
    os << header("UTF-8")
       << stylesheet(xslt_path("ALPS.xsl"))
       << start_tag("SIMULATION")
       << xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
       << attribute("xsi:noNamespaceSchemaLocation", "http://xml.comp-phys.org/2002/10/ALPS.xsd");
    os << params;
    if (write_xml) {
      if (obs.size() == 1) {
        scheduler::write_xml_averages(os, obs[0]);
      } else {
        for (std::size_t i = 0; i < obs.size(); ++i) scheduler::write_xml_averages(os, obs[i], i+1);
      }
    }
    for (std::size_t i = 0; i < info.size(); ++i) os << info[i];
    os << end_tag("SIMULATION");
    if (make_backup && exists(file_bak)) remove(file_bak);
  }
};


// Reads the parameters and the clones of a task. The results are not read
// from the XML file: they are in the HDF5 files of the task.
class simulation_xml_handler : public CompositeXMLHandler {
public:
  simulation_xml_handler(Parameters& params, std::deque<clone_info>& clones) :
    CompositeXMLHandler("SIMULATION"), params_handler_(params), os_handler_("AVERAGES"),
    clones_(clones), clone_handler_(clone_) {
    add_handler(params_handler_);
    add_handler(os_handler_);
    add_handler(clone_handler_);
  }

protected:
  void start_child(const std::string& name, alps::XMLAttributes const& /* attributes */,
    xml::tag_type type) {
    if (type == alps::xml::element && name == "MCRUN")
      clone_ = clone_info();
  }

  void end_child(const std::string& name, xml::tag_type type) {
    if (type == alps::xml::element && name == "MCRUN")
      clones_.push_back(clone_);
  }

private:
  ParametersXMLHandler params_handler_;
  DummyXMLHandler os_handler_;
  std::deque<clone_info>& clones_;
  clone_info clone_;
  clone_info_xml_handler clone_handler_;
};

class simulation_parameters_xml_handler : public CompositeXMLHandler {
public:
  simulation_parameters_xml_handler(Parameters& params) :
    CompositeXMLHandler("SIMULATION"), params_handler_(params), os_handler_("AVERAGES"),
    clone_handler_("MCRUN") {
    add_handler(params_handler_);
    add_handler(os_handler_);
    add_handler(clone_handler_);
  }

private:
  ParametersXMLHandler params_handler_;
  DummyXMLHandler os_handler_;
  DummyXMLHandler clone_handler_;
};

} // end namespace alps

#endif // PARAPACK_SIMULATION_P_H
