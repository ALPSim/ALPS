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
#ifndef HYB_SEG_HPP
#define HYB_SEG_HPP
#include<fstream>

//the struct 'segment' has a start and an end time and a comparison function.
//segment describes a pair of operators in an orbital: t_start is the c^\dagger, t_end is the c operator.
struct segment
{
  segment(){};
  segment(double t_start, double t_end){t_start_=t_start; t_end_=t_end;}
  double t_start_, t_end_;
  
  bool operator<(const segment &t2)const { return t_start_<t2.t_start_; }
  bool operator<(const double &t2) const { return t_start_<t2; }
};

std::ostream &operator<<(std::ostream &os, const segment &s);

#endif
