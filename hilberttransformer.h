/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_DMFT_HILBERTTRANSFORMER_H
#define ALPS_DMFT_HILBERTTRANSFORMER_H

#include "solver.h"
#include <functional>


class HilbertTransformer : public std::binary_function<const vector_type&,double,vector_type>
{
public:
  virtual result_type operator()(first_argument_type Sigma, second_argument_type mu) =0; 
  virtual ~HilbertTransformer() {}
};

class SemicircleHilbertTransformer : public HilbertTransformer {
public:
  SemicircleHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  result_type operator()(first_argument_type Sigma, second_argument_type mu)
  {
   // to implement
  return vector_type();
  
  }
private:
  double bw_; // or band width
};

class SquarelatticeHilbertTransformer : public HilbertTransformer {
public:
  SquarelatticeHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  result_type operator()(first_argument_type Sigma, second_argument_type mu)
  {
    //to implement
  
  return vector_type();
  }
private:
  double bw_; // or band width
};

// ...

class GeneralHilbertTransformer : public HilbertTransformer {
public:
  GeneralHilbertTransformer(const vector_type& D) : D_(D) {}
  result_type operator()(first_argument_type Sigma, second_argument_type mu)
  {
    //to implement
  
  return vector_type();
  }
private:
  vector_type D_; // or band width
};


#endif
