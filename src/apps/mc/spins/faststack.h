/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
