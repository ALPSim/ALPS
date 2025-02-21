/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef ALPS_MPSPSCAN_SCHEDULER_HPP
#define ALPS_MPSPSCAN_SCHEDULER_HPP

#include "dmrg/utils/time_stopper.h"
#include "libpscan/options.hpp"
#include <boost/filesystem.hpp>

enum TaskStatusFlag {
    TaskNotStarted, TaskRunning, TaskHalted, TaskFinished
};

struct TaskDescriptor {
    TaskStatusFlag status;
    boost::filesystem::path in;
    boost::filesystem::path out;
};

class Scheduler {
public:
    Scheduler(const Options&);
    void run();
    
private:
    void parse_job_file(const boost::filesystem::path&);
    void checkpoint_status() const;
    
    boost::filesystem::path outfilepath;
    boost::filesystem::path infilepath;
    
    time_stopper stop_callback;
    bool write_xml;
    std::vector<TaskDescriptor> tasks;
};

#endif
