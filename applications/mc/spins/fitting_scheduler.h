/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_FITTINGSCHEDULERS_H
#define ALPS_APPLICATIONS_FITTINGSCHEDULERS_H

#include <alps/utility/copyright.hpp>
#include <alps/scheduler.h>
#include <alps/scheduler/factory.h>
#include <alps/scheduler/options.h>
#include <alps/scheduler/types.h>
#include <alps/scheduler/signal.hpp>
#include <boost/lexical_cast.hpp>

#include <alps/plot.h>

#include <cstdio>
#include <iostream>
#include <fstream>

#include "fitter.h"
#include <alps/scheduler/options.h>

namespace alps {
namespace scheduler {
//<<<<<<< fitting_scheduler.h

template <class PARENT>
class FittingScheduler : public PARENT
{
public:
  /**
   * Constructor. Read parameters and determine the type of fitter to be used.
   */
  FittingScheduler(const NoJobfileOptions& opt, const Factory& p, 
               const char* sim_file_name) 
      : PARENT(opt,p)
  { 
    std::cerr << "\nFittingScheduler is being constructed\n"; 
    checkpoint_filename = boost::filesystem::absolute(
               std::string(sim_file_name) + ".checkpoint");

    alps::Parameters sim_param;
    {
      std::ifstream in(sim_file_name);
      in >> sim_param;
    }

    if (sim_param.defined("FITTER")) {
      if (sim_param["FITTER"] == "Dummy")
        theFitter = new DummyFitter(sim_file_name);
      else if (sim_param["FITTER"] == "EstGradient")
        theFitter = new EstGradFitter(sim_file_name);
      else {
        std::cerr << "Invalid parameter for the type of the fitter:\n"
                  << sim_param["FITTER"] << " is not a valid option\n";
        boost::throw_exception(std::runtime_error("invalid fitter specified"));
      }
    } else
      theFitter = new EstGradFitter(sim_file_name);

    namebase = theFitter->get_filenamebase();
    fit_observable = theFitter->get_fit_observable();
    exp_results_set = theFitter->get_exp_results_plot_set();
    exp_error_set = theFitter->get_exp_error_plot_set();
    exp_percent_set = theFitter->get_exp_percent_plot_set();

    loaded_from_checkpoint = false;
    if (boost::filesystem::exists(checkpoint_filename)) {
      // read from checkpoint file
      std::ifstream in(checkpoint_filename.filename().string().c_str());
      theFitter->initialize(in);
      loaded_from_checkpoint = true;
    } else
      theFitter->initialize();
    PARENT::make_summary = true;

    if (sim_param.defined("TIME_LIMIT"))
      PARENT::time_limit = (double)sim_param["TIME_LIMIT"];

    // options for the plots of the intermediary fitting results
    make_results_plots = (sim_param.defined("MAKE_RESULTS_PLOTS")
                        && (sim_param["MAKE_RESULTS_PLOTS"] == "yes"));
    if (make_results_plots) {
      exp_results_set = theFitter->get_exp_results_plot_set();
      std::cout << "Result plot will be made.\n";
    }      
    show_results_plots = (make_results_plots 
                        && sim_param.defined("SHOW_RESULTS_PLOTS")
                        && (sim_param["SHOW_RESULTS_PLOTS"] == "yes"));

    make_residui_plots = (sim_param.defined("MAKE_RESIDUI_PLOTS")
                        && (sim_param["MAKE_RESIDUI_PLOTS"] == "yes"));
    if (make_residui_plots) {
      exp_error_set = theFitter->get_exp_error_plot_set();
      std::cout << "Residui plot will be made.\n";
    }
    show_residui_plots = make_residui_plots
                        && (sim_param.defined("MAKE_RESIDUI_PLOTS")
                        && (sim_param["MAKE_RESIDUI_PLOTS"] == "yes"));

    make_percent_plots = (sim_param.defined("MAKE_PERCENT_PLOTS")
                        && (sim_param["MAKE_PERCENT_PLOTS"] == "yes"));
    if (make_percent_plots) {
      exp_percent_set = theFitter->get_exp_percent_plot_set();
      std::cout << "Percent plots will be made.\n";
    }
    show_percent_plots = make_percent_plots
                        && (sim_param.defined("SHOW_PERCENT_PLOTS")
                        && (sim_param["SHOW_PERCENT_PLOTS"] == "yes"));

    // options involving all plot types
    if (sim_param.defined("MAKE_PLOTS")
               && (sim_param["MAKE_PLOTS"] == "yes")) {
      make_results_plots = true;
      make_residui_plots = true;
      make_percent_plots = true;
      std::cout << "All plot types will be made.\n";
    }
    if (sim_param.defined("SHOW_PLOTS")
               && (sim_param["SHOW_PLOTS"] == "yes")) {
      show_results_plots = true;
      show_residui_plots = true;
      show_percent_plots = true;
    }

    // options for the file handling
    delete_runfiles = (sim_param.defined("DELETE_RUNFILES")
                        && (sim_param["DELETE_RUNFILES"] == "yes"));
    if (delete_runfiles)
      std::cout << "run files will be deleted\n";
    else
      std::cout << "run files will not be deleted\n";

    delete_outfiles = (sim_param.defined("DELETE_OUTFILES") 
                        && (sim_param["DELETE_OUTFILES"] == "yes"));
    if (delete_outfiles)
      std::cout << "output files will be deleted\n";
    else
      std::cout << "output files will not be deleted\n";

    delete_infiles = (sim_param.defined("DELETE_INFILES")
                        && (sim_param["DELETE_INFILES"] == "yes"));
    if (delete_infiles)
      std::cout << "input files will be deleted\n";
    else
      std::cout << "input files will not be deleted\n";

    delete_some_files= ( delete_runfiles || delete_outfiles || delete_infiles );
  }

  /** 
   * The main loop, alternating between the generation of new job files and 
   * passing the control to the superclass-function run.
   */
  int run() {
    std::string curFileName;
    double seconds_left;
    while (!(theFitter->is_done(&seconds_left))) {
      std::cout << "starting next run, " << PARENT::time_limit
                << " seconds left\n";
      if (loaded_from_checkpoint)
        loaded_from_checkpoint = false;
      else
        theFitter->make_new_jobfile();
      curFileName = theFitter->get_current_filename();
      std::cerr << "new job file name is " << curFileName << "\n";
      theScheduler->set_new_jobfile(curFileName);

      if (seconds_left > 0.)
        PARENT::time_limit = seconds_left;
      std::cout << "starting new MC job\n";
      int res = PARENT::run();

      if (res<0) {
        std::cerr << "computation interrupted, returning ... \n";
        end_fitting();
        return 0;
      }
      theFitter->add_sim_result(theScheduler->getSummary());
      theFitter->print_results();
 
      if (make_results_plots) 
        make_results_plot(theFitter->get_step_count());
      if (make_residui_plots)
        make_residual_plot(theFitter->get_step_count());
      if (make_percent_plots)
        make_percent_plot(theFitter->get_step_count());

      if (delete_some_files)
        delete_files();
    }
    end_fitting();
  return 0;
  }

protected:
  AbstractFitter* theFitter;
  boost::filesystem::path checkpoint_filename;
  bool loaded_from_checkpoint;

  bool make_results_plots, show_results_plots;
  bool make_residui_plots, show_residui_plots;
  bool make_percent_plots, show_percent_plots;

  bool delete_runfiles, delete_outfiles, delete_infiles, delete_some_files;
  std::string fit_observable;
  std::string namebase;

  alps::plot::Set<double> exp_results_set, exp_error_set, exp_percent_set;

  /**
   * Make a checkpoint, mainly by calling the function of the fitter.
   */
  virtual void checkpoint() {
    bool make_backup=boost::filesystem::exists(checkpoint_filename);
    boost::filesystem::path fname = checkpoint_filename;
    boost::filesystem::path dir = checkpoint_filename.parent_path();
                                                                                
    if (make_backup)
      fname=dir/(fname.filename().string()+".bak");
                                                                                
    oxstream out(fname);
    out << header("UTF-8") << stylesheet(xslt_path("ALPS.xsl"));
    theFitter->write_to_xml(out);
                                                                                
    if (make_backup) {
      boost::filesystem::remove(checkpoint_filename);
      boost::filesystem::rename(fname,checkpoint_filename);
    }
    // checkpoint the tasks, as it is done in the parent class
    PARENT::checkpoint();
  }

  /** 
   * Make a results plot.
   */
  void make_results_plot(int step) const {
    alps::plot::Set<double> sim_results_set;
    sim_results_set = theFitter->get_sim_results_plot_set();

    char title[50];
    snprintf(title, sizeof(title), "Fitting after %i steps", step);
    alps::plot::Plot<double> myPlot(title);
    myPlot << exp_results_set;
    myPlot << sim_results_set;
    
    // label for x and y axis
    myPlot.set_labels("Temperature [K]", fit_observable);
    myPlot.show_legend(true);

    char plot_xml_name[128], plot_xmgr_name[128];
    snprintf(plot_xml_name, sizeof(plot_xml_name), "%s_results_plot%i.xml",namebase.c_str(),step);
    snprintf(plot_xmgr_name, sizeof(plot_xmgr_name), "%s_results_plot%i.xmgr", namebase.c_str(),step);
    make_and_show(myPlot, plot_xml_name, plot_xmgr_name, show_results_plots);
  }

  /**
   * Make a residual plot of the current step.
   */
  void make_residual_plot(int step) const {
    alps::plot::Set<double> residui_set;
    residui_set = theFitter->get_res_abs_plot_set();

    char title[50];
    snprintf(title, sizeof(title), "Residui after %i Steps and Measurement Error", step);
    alps::plot::Plot<double> myPlot(title);
    myPlot << exp_error_set;
    myPlot << residui_set;

    // label the axis
    myPlot.set_labels("Temperature [K]", fit_observable);
    myPlot.show_legend(true);

    char plot_xml_name[128], plot_xmgr_name[128];
    snprintf(plot_xml_name, sizeof(plot_xml_name), "%s_residual_plot%i.xml", namebase.c_str(), step);
    snprintf(plot_xmgr_name, sizeof(plot_xmgr_name), "%s_residual_plot%i.xmgr", namebase.c_str(), step);
    make_and_show(myPlot, plot_xml_name, plot_xmgr_name, show_residui_plots);
  }

  /**
   * Make a percent plot.
   */
  void make_percent_plot(int step) const {
    alps::plot::Set<double> sim_percent = theFitter->get_res_percent_plot_set();
    
    char title[64];
    snprintf(title, sizeof(title), "Residuum Size and Measurement Error after %i steps", step);
    alps::plot::Plot<double> myPlot(title);
    myPlot << exp_percent_set;
    myPlot << sim_percent;
    
    char plot_xml_name[128], plot_xmgr_name[128];
    snprintf(plot_xml_name, sizeof(plot_xml_name), "%s_percent_plot%i.xml", namebase.c_str(), step);
    snprintf(plot_xmgr_name, sizeof(plot_xmgr_name), "%s_percent_plot%i.xmgr", namebase.c_str(), step);
    
    myPlot.set_labels("Temperature [K]", "Percents of Measurement value");
    myPlot.show_legend(true);
    make_and_show(myPlot, plot_xml_name, plot_xmgr_name, show_percent_plots);
  }

  /**
   * Function to convert xml to xmgr files and to possibly display plots using 
   * grace.
   */
  void make_and_show(const alps::plot::Plot<double>& plot,
               const char xml_name[], const char xmgr_name[], bool show) const {
    alps::oxstream out(xml_name);
    out << plot;
                                                                                
    char make_plot[128];
    snprintf(make_plot, sizeof(make_plot), "plot2xmgr %s > %s", xml_name, xmgr_name);
    system(make_plot);
                                                                                
    if (show) {
      char show_plot[64];
      snprintf(show_plot, sizeof(show_plot), "xmgrace %s &", xmgr_name);
      system(show_plot);
    }
  }
                                                                                
  /** 
   * Gracefully terminate the fitting by writing down all the results.
   */
  virtual void end_fitting() {
    theFitter->print_results();
    delete theFitter;
  }
     
  /** 
    * Delete files that are not needed any more by the fitting classes.
    */
  virtual void delete_files() {
    int step_count = theFitter->get_step_count();
    char cmd[256];

    if (delete_infiles) {
      snprintf(cmd, sizeof(cmd), "rm %s.step%i.in.xml %s.step%i.task*.in.xml",
               namebase.c_str(),step_count, namebase.c_str(), step_count);
      system(cmd);
    }

    if (delete_runfiles) {
      // delete the runfiles of the previous step, as the ones of the current
      snprintf(cmd, sizeof(cmd), "rm %s.step%i.task*.out.run*", namebase.c_str(),step_count-1);
      system(cmd);
    }

    if (delete_outfiles) {
      // again, delete the files of the previous step
      snprintf(cmd, sizeof(cmd), "rm %s.step%i.out.xml %s.step%i.task*.out.xml",
               namebase.c_str(),step_count-1, namebase.c_str(), step_count-1);
      system(cmd);
    }
  }
};
                                                                                

int start_fitting(int argc, char** argv, const Factory& p);

//=======
//  int sim_start(int argc, char** argv, const Factory& p);
//  void end_fitting(bool isManager);
// the fitter on this node as global variable
//  extern AbstractFitter* theFitter;
//>>>>>>> 1.2
} // namespace scheduler
} // namespace alps

#endif // ALPS_APPLICATIONS_FITTINGSCHEDULERS_H
