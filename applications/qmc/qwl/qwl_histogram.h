/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004 by Stefan Wessel <wessel@comp-phys.org>
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

#ifndef ALPS_QWL_HISTOGRAM_H
#define ALPS_QWL_HISTOGRAM_H

#include "alps/scheduler/montecarlo.h"
#include <vector>
#include <algorithm>
#include <valarray>

template<typename T=double> 
class histogram {
 public:
  histogram() : data(), left_(0), right_(0) {}
  void resize(unsigned int newsize,unsigned int newleft=0) {
    data.resize(newsize);
    left_=newleft;
    right_=left_+size()-1;
    fill(0);
  }
  void fill(T x) {
    std::fill(data.begin(), data.end(), x);
  }
  void subtract() {
    T x=data[0];
    for (unsigned int i=0;i<data.size();++i)
      data[i]-=x;
  }
  std::valarray<T> getvalarray(unsigned int min,unsigned int max) const {
    std::valarray<T> dval;
    dval.resize(max-min+1);
    int count=0;
    for (int i=min-left();i<=max-left();++i) {
      dval[count]=data[i];
      ++count;
    }
    return dval;
  }
  T& operator[](unsigned int i) {
    return data[i-left()];
  }
  unsigned int left() const {
    return left_;
  }  
  unsigned int right() const {
    return right_;
  }  
  unsigned int size() const {
    return data.size();
  }
  double flatness() const {  
    double av=data[0];
    for (int i=1;i<size();++i)
      av+=data[i];
    av/=size();
    double diff=fabs(data[0]-av);
    for (int i=1;i<size();++i)
      if (diff<fabs(data[i]-av))  
        diff=fabs(data[i]-av);
    return ( av>0 ? diff/av : 0);
  }  
  T min() const {  
   T min=data[0];
    for (int i=1;i<size();++i)
      if (data[i]<min)
         min=data[i];
    return min;
  }  
  void save(alps::ODump& dump) const {
    dump << data << left_;
  }
  void load(alps::IDump& dump) {
    dump >> data >> left_;
    right_=left_+size()-1;
  }
 private:
  std::vector<T> data;
  unsigned int left_;
  unsigned int right_;
};

#endif


