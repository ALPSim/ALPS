/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
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

#include <alps/ngs.hpp>
#include <fstream>
#include "../green_function.h"
#ifndef HYB_FUN_HPP
#define HYB_FUN_HPP

//container for the hybridization function Delta(\tau) (historically: F(\tau) )
class hybfun : public green_function<double>{
  public:
    //constructor
  hybfun(const alps::params &p);
  hybfun(const hybfun &rhs):green_function<double>(rhs),beta_(rhs.beta_) {
//    operator=(rhs);
//    for (int nf=0;nf<nflavor();nf++)
//      for (int nt=0;nt<ntime();nt++)
//        operator()(nt,nf)=rhs(nt,nf);
//    std::cerr << rhs.ntime() << " " << ntime() << " " << rhs.nflavor() << " " << nflavor() << std::endl;
    hybridization_function_sanity_check();
  }
  const hybfun &operator=(const hybfun &rhs) {
    beta_ = rhs.beta_;
//    green_function<double>::operator=(rhs);
//    std::cerr << rhs.ntime() << " " << ntime() << " " << rhs.nflavor() << " " << nflavor() << std::endl;
    for (int nf=0;nf<nflavor();nf++)
      for (int nt=0;nt<ntime();nt++)
        operator()(nt,nf)=rhs(nt,nf);
    hybridization_function_sanity_check();
    return *this;
  }
  
  ~hybfun() {
//    std::cerr << "Deleting hybfunc\n";
  }
  
  double interpolate(double time, int orbital) const;

  friend std::ostream &operator<<(std::ostream &os, const hybfun &hyb);
private:
  void read_hybridization_function(const alps::params &p);
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

#endif
