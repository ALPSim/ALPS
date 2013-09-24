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

#include<fstream>

//this keeps track of values for the bare green's function at times
//tau_A-tau_i, where i is a vertex and A is a measuring point for the imag
//time greens function.

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
    if(new_nop<=(unsigned int)nop_){ //down is easy
      nop_=new_nop;
      return;
    } else if(new_nop<= (unsigned int)memory_size_){ //up is easy as long as we don't have to allocate new memory
      nop_=new_nop;
    } else{ //get new memory
      double *new_values_=new double[nt_*new_nop];
      memcpy(new_values_, values_, sizeof(double)*nop_*nt_);
      delete[] values_;        //free memory
      values_=new_values_;        //let the matrix point to the new memory location.
      nop_=new_nop;
      memory_size_=new_nop;
    }
  }
  
  inline void swap_vertices(unsigned int p, unsigned int q)
  {
    fortran_int_t inc=1;
    FORTRAN_ID(dswap)(&nt_, values_+p*nt_, &inc, values_+q*nt_, &inc);
  }
  
  inline double* values(){return values_;}
  inline int memory_size(){return memory_size_;}

private:
  fortran_int_t memory_size_;
  fortran_int_t nop_;
  fortran_int_t nt_;
  double *values_;
};
