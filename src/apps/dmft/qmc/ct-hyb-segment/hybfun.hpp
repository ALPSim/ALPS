/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <fstream>
#include <alps/params.hpp>
#include <alps/gf/gf.hpp>
#include <alps/gf/tail.hpp>
#include <alps/gf/fourier.hpp>
#include "./green_function.h"

#ifndef HYB_FUN_HPP
#define HYB_FUN_HPP

// Internal to ALPS::cthyb_segment: keeps its classes apart from ALPS::cthyb's.
namespace cthyb_segment {


//container for the hybridization function Delta(\tau) (historically: F(\tau) )
class hybfun : public green_function<double>{
  public:
    //constructor
  hybfun(const alps::params &p);
  double interpolate(double time, int orbital) const;

  friend std::ostream &operator<<(std::ostream &os, const hybfun &hyb);
private:
  void read_hybridization_function(const alps::params &p);
  void read_hybridization_from_h5gf(const alps::params &p);
  void hybridization_function_sanity_check(void);
  double beta_;
};

std::ostream &operator<<(std::ostream &os, const hybfun &hyb);

#ifndef COLORS
#define COLORS
#define cblack "\033[22;30m"
#define cred "\033[22;31m"
#define cgreen "\033[22;32m"
#define cbrown "\033[22;33m"
#define cblue "\033[22;34m"
#define cmagenta "\033[22;35m"
#define ccyan "\033[22;36m"
#define cgray "\033[22;37m"
#define cdgray "\033[01;30m"
#define clred "\033[01;31m"
#define clgreen "\033[01;32m"
#define clyellow "\033[01;33m"
#define clblue "\033[01;34m"
#define clmagenta "\033[01;35m"
#define clcyan "\033[01;36m"
#define cwhite "\033[01;37m"
#endif


}  // namespace cthyb_segment
#endif
