/*****************************************************************************
*
* Copyright (C) 2005 - 2015 Emanuel Gull <egull@umich.edu>
*
*****************************************************************************/
#pragma once
#include <complex>
#include <vector>
#include <iostream>

///enum for spin up and spin down
enum  {up=0, down=1} ;

struct hifreq_moments{
  std::vector<std::vector<std::vector<double> > >c0;
  std::vector<std::vector<std::vector<double> > >c1;
  std::vector<std::vector<std::vector<double> > >c2;
  std::vector<std::vector<std::vector<double> > >c3;
};

namespace alps{class Parameters;}
