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
