/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#include<fstream>

//this keeps track of values for the bare green's function at times
//tau_A-tau_i, where i is a vertex and A is a measuring point for the imag
//time greens function.

extern "C" void dswap_(unsigned int* n, double *a, int *inca, double *b, int *incb);
double green0_spline(const double, const unsigned int,const unsigned int,const unsigned int);

class green_matrix
{
public:
  green_matrix(const unsigned int n_tau, const unsigned int noperators)
  {
    values_=new double[(n_tau+1)*noperators];
    nop_=noperators;
    memory_size_=noperators;
    nt_=n_tau+1; //plus one because we want to take into account the last point
  }

  ~green_matrix()
  {
    delete[] values_;
  }
	
  const green_matrix& operator=(const green_matrix &g)
  {
    memory_size_=g.memory_size_;
    nop_=g.nop_;
    nt_=g.nt_;
    memcpy(values_, g.values_, memory_size_*nt_*sizeof(double));
    return *this;
  }
  
  green_matrix(const green_matrix &g)
  {
    values_=new double[g.memory_size_*g.nt_];
    operator=(g);
  }
  
  inline double &operator()(const int tau, const int op)
  {
    return *(values_+op*nt_+tau);
  }
  
  inline const double &operator()(const int tau, const int op) const 
  {
    return *(values_+(op*nt_+tau));
  }
  
  void resize(const unsigned int new_nop)
  {
    if(new_nop<=nop_){ //down is easy
      nop_=new_nop;
      return;
    } else if(new_nop<= memory_size_){ //up is easy as long as we don't have to allocate new memory
      nop_=new_nop;
    } else{ //get new memory
      double *new_values_=new double[nt_*new_nop];
      memcpy(new_values_, values_, sizeof(double)*nop_*nt_);
      delete[] values_;	//free memory
      values_=new_values_;	//let the matrix point to the new memory location.
      nop_=new_nop;
      memory_size_=new_nop;
    }
  }
  
  inline void swap_vertices(unsigned int p, unsigned int q)
  {
    int inc=1;
    dswap_(&nt_, values_+p*nt_, &inc, values_+q*nt_, &inc);
  }
  
  inline double* values(){return values_;}
  inline int memory_size(){return memory_size_;}

private:
  unsigned int memory_size_;
  unsigned int nop_;
  unsigned int nt_;
  double *values_;
};
