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

#ifndef DMFT_QMC_WEAK_COUPLING_RESIZABLE_MATRIX_H
#define DMFT_QMC_WEAK_COUPLING_RESIZABLE_MATRIX_H

#include<fstream>
#include "operator.hpp"
#include <boost/numeric/bindings/blas.hpp>
//this implements all functinons necessary for an inverse 'M' - matrix.

double inner_prod(const double* v1, const double *v2, const int size);
void scale(const double alpha, double *v, const fortran_int_t size);


inline double vector_max(const double *v, const fortran_int_t noperators)
{
  fortran_int_t inc=1;
  fortran_int_t max_index=FORTRAN_ID(idamax)(&noperators, v ,&inc);
  return fabs(v[max_index-1]); //fortran indexing
}


class resizeable_matrix
{
public:

  resizeable_matrix(const unsigned int ms=1)
  {
    values_=new double[ms*ms];
    size_=0;
    memory_size_=ms;
  }

  ~resizeable_matrix()
  {
    delete[] values_;
  }
  
  const resizeable_matrix& operator=(const resizeable_matrix &M)
  {
    memory_size_=M.memory_size_;
    size_=M.size_;
    max_size_=M.max_size_;
    memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(double));
    creators_=M.creators_;
    annihilators_=M.annihilators_;
    alpha_=M.alpha_;
    return *this;
  }
  
  resizeable_matrix(const resizeable_matrix &M)
  {
    values_=new double[M.memory_size_*M.memory_size_];
    operator=(M);
  }
  
  inline double &operator()(const unsigned int i, const unsigned int j){return *(values_+(i*memory_size_+j));}
  inline const double &operator()(const unsigned int i, const unsigned int j) const {return *(values_+(i*memory_size_+j));}

  void resize(const unsigned int new_size)
  {
    if(new_size<=size_){ //down is easy
      if((int)new_size < (int)(size_-5) && (int) new_size > 10){
        double *new_values=new double[new_size*new_size];
        for(unsigned int i=0;i<size_;++i){
          memcpy(new_values+i*new_size, values_+i*memory_size_, sizeof(double)*new_size); //for each row: copy the entire row.
        }
        delete[] values_;       //free memory
        values_=new_values;    //let the matrix point to the new memory location.
        size_=new_size;
        memory_size_=new_size;
      }
      size_=new_size;
      return;
    } else if(new_size<= memory_size_){ //up is easy as long as we don't have to allocate new memory
      size_=new_size;
    } else{ //get new memory
      double *new_values=new double[new_size*new_size];
      for(unsigned int i=0;i<size_;++i){
        memcpy(new_values+i*new_size, values_+i*memory_size_, sizeof(double)*size_); //for each row: copy the entire row.
      }
      delete[] values_;    //free memory
      values_=new_values;    //let the matrix point to the new memory location.
      size_=new_size;
      memory_size_=new_size;
    }
  }
  
  inline unsigned int size()const{return size_;}
  
  inline void add_outer_prod(const double *v1, const double *v2)
  {
  //smart people call this 'rank one update'.
    double alpha=1.;
    fortran_int_t inc=1;
    FORTRAN_ID(dger)(&size_, &size_, &alpha, v2, &inc, v1, &inc, values_, &memory_size_); 
  }

  inline void right_multiply(const double *v1, double *v2)
  { //perform v2[i]=M[ij]v1[j]
    //call the BLAS routine for matrix vector multiplication:
    char trans='T';
    double alpha=1., beta=0.;    //no need to multiply a constant or add a vector
    int inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1, &inc, &beta, v2, &inc);
  }


  inline void matrix_right_multiply(const double *M1, double *M2, const fortran_int_t columns)

  { //perform M2(i,j)=sum_k M(i,k)*M1(k,j)
    //call the BLAS routine for matrix matrix multiplication:     
    char transa='N';
    char transb='N';
    double alpha=1.;
    double beta=0.;
    FORTRAN_ID(dgemm)(&transa, &transb, &columns, &size_, &size_, &alpha, M1, &columns, values_, &memory_size_, &beta, M2, &columns);
  }
  
  void left_multiply(const double *v1, double *v2)
  { //perform v2[i]=v1[j]M[ji]
    //call the BLAS routine for matrix vector multiplication:
    char trans='N';
    double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
    int inc=1;
    FORTRAN_ID(dgemv)(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1, &inc, &beta, v2, &inc);
  }
  
  bool check()
  {
    bool pass=true;
    for(unsigned int i=0;i<size_;++i){
      for(unsigned int j=0;j<size_;++j){
        if(operator()(i,j)>2){
          std::cout<<"M is really large: "<<i<<" "<<j<<" "<<operator()(i,j)<<" "<<size_<<std::endl;
          pass=false;
        }
      }
    }
    return pass;
  }
  
  double max() const
  {

    double* rowmax = new double[size_];

    fortran_int_t rowmax_index, max_index, inc=1;
    for(unsigned int i=0;i<size_;++i){
      rowmax_index=FORTRAN_ID(idamax)(&size_, values_+i*memory_size_,&inc);
      rowmax[i]=*(values_+i*memory_size_+rowmax_index-1); //fortran convention: start counting from one
    }
    max_index=FORTRAN_ID(idamax)(&size_, rowmax ,&inc);

    double m=fabs(rowmax[max_index-1]);
    delete[] rowmax;
    return m;

  }
  
  // return functions
  void clear(){ size_=0;}
  inline unsigned int memory_size() const {return memory_size_;}
  
  std::vector<creator> &creators(){ return creators_;}
  const std::vector<creator> &creators() const{ return creators_;}
  std::vector<annihilator> &annihilators(){ return annihilators_;}
  const std::vector<annihilator> &annihilators()const{ return annihilators_;}
  std::vector<double> &alpha(){ return alpha_;}
  const std::vector<double> &alpha() const{ return alpha_;}
  
  
private:
  std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix
  std::vector<annihilator> annihilators_; //an array of to annihilation operators c corresponding to the column of the matrix
  std::vector<double> alpha_;             //an array of doubles corresponding to the alphas of Rubtsov for the c, cdaggers at the same index.
  
  fortran_int_t memory_size_; //size of matrix for which we have memory allocated
  fortran_int_t size_; //current size of matrix
  fortran_int_t max_size_; //max size of matrix
  double *values_; //where the actual values are stored
};

std::ostream & operator<<(std::ostream &os, const resizeable_matrix &M);

#endif
