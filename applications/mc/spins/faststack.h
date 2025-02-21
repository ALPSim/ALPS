/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>
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

#ifndef ALPS_APPLICATIONS_MC_SPIN_FASTSTACK_H_
#define ALPS_APPLICATIONS_MC_SPIN_FASTSTACK_H_

template <class T> class fast_stack {
public:
  fast_stack(std::size_t max_size) 
    : stack_(new T[max_size]), ptr_(stack_-1), start_(stack_-1) 
    {}
  ~fast_stack() { delete[] stack_;}
  T& top() { return *ptr_;}
  void pop() { --ptr_;}
  void push(T x) {*(++ptr_) = x;}
  bool empty() { return ptr_==start_;}
private:
  T *stack_, *ptr_, *start_;
};

#endif
