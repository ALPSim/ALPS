#include <functional>
#include <solver.h>


class HilbertTransformer : public std::binary_function<const vector_type&,double,vector_type>
{
public:
  virtual result_type operator(first_argument_type Sigma, second_argument_type mu) =0; 
  virtual ~HilbertTransformer() {}
}

class SemicircleHilbertTransformer : public HilbertTransformer {
public:
  SemicircleHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  result_type operator(first_argument_type Sigma, second_argument_type mu)
  {
    to implement
  
  
  }
private:
  double bw_; // or band width
}

class SquarelatticeHilbertTransformer : public HilbertTransformer {
public:
  SquarelatticeHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  result_type operator(first_argument_type Sigma, second_argument_type mu)
  {
    to implement
  
  
  }
private:
  double bw_; // or band width
}

..

class GeneralHilbertTransformer : public HilbertTransformer {
public:
  GeneralHilbertTransformer(const vector_type& D) : D_(D) {}
  result_type operator(first_argument_type Sigma, second_argument_type mu)
  {
    to implement
  
  
  }
private:
  vector_type D_; // or band width
}

